
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "vcf.h"

#include "PhaseImpMissing.hpp"
#include "FileUtils.hpp"
#include "format.hpp"
#include "timestamp.hpp"

using namespace std;
using namespace Osprey;

extern vector< vector<double> > phaseImpCore(const vector< vector< pair<double, int> > >& ibsMatrix,
                                             const vector< vector<double> >& dipCNPs,
                                             const int nIterations,
                                             const int debug);

static const double MAX_GENOTYPE_QUALITY = 99.0;
static const double MIN_CN_LIKELIHOOD = -1000.0;

static const double weight0err = 0.5;
static const double weightTotalLen = 0.25;

namespace Osprey {

    PhaseImpMissing::PhaseImpMissing() {
        mDebug = 0;
        mVerbose = 0;
        mIterations = 0;
        mOutputHeader = NULL;
    }

    PhaseImpMissing::~PhaseImpMissing() {
    }

    void PhaseImpMissing::setDebug(int value) {
        mDebug = value;
    }

    void PhaseImpMissing::setVerbose(int value) {
        mVerbose = value;
    }

    void PhaseImpMissing::setIterations(int value) {
        mIterations = value;
    }

    static void addHeaderLine(bcf_hdr_t* header, const char* text) {
        int len = 0;
        bcf_hrec_t* hrec = bcf_hdr_parse_line(header, text, &len);
        if (hrec == NULL) {
            throw std::runtime_error(format("Internal error: Failed to parse header line: %s", text));
        }
        if (bcf_hdr_add_hrec(header, hrec) < 0) {
            throw std::runtime_error(format("Internal error: Failed to add header line: %s", text));
        }
    }

    void PhaseImpMissing::updateVCFHeader(VCFReader& vcfReader, VCFWriter& vcfWriter) {
        bcf_hdr_t* header = bcf_hdr_dup(vcfReader.getHeader());
        addHeaderLine(header, "##FORMAT=<ID=PCN,Number=1,Type=String,Description=\"Phased copy number\">");
        addHeaderLine(header, "##FORMAT=<ID=PCNF,Number=1,Type=String,Description=\"Phased fractional copy number\">");
        addHeaderLine(header, "##FORMAT=<ID=PCNQ,Number=1,Type=String,Description=\"Phased copy number quality\">");
        addHeaderLine(header, "##FORMAT=<ID=PCNL,Number=1,Type=String,Description=\"Phased copy number likelihoods\">");
        addHeaderLine(header, "##FORMAT=<ID=PST,Number=1,Type=String,Description=\"Phasing status (PHASED/IMPUTED)\">");

        int nImputed = 0;
        int nSamples = mSampleIds.size();
        mSampleStatus.resize(nSamples, string("."));
        const vector<string> inputSamples = vcfReader.getSampleIds();
        const set<string> inputSampleSet(inputSamples.begin(), inputSamples.end());
        for (int i = 0; i < nSamples; i++) {
            const string& sample = mSampleIds[i];
            if (inputSampleSet.find(sample) == inputSampleSet.end()) {
                mSampleStatus[i] = "I";
                bcf_hdr_add_sample(header, sample.c_str());
                nImputed++;
            } else {
                mSampleStatus[i] = "P";
            }
        }
        int status = bcf_hdr_sync(header);
        if (status != 0) {
            throw std::runtime_error(format("Error in bcf_hdr_sync: status code %d", status));
        }

        vcfWriter.writeHeader(header);
        mOutputHeader = header;

        if (mVerbose > 0) {
            cout << timestamp() << " Found " << nSamples << " samples in input"
                 << ", " << nImputed << " will be imputed" << endl;
        }
    }

    void PhaseImpMissing::readIBSMatrix(std::string inputFile) {

        // Clear the current matrix state
        mSampleIds.clear();
        mSampleIndexMap.clear();
        mIBSMatrix.clear();

        // First pass over input file: read sample IDs
        // TBD: Not tab-delimited?
        FileUtils::AutoGzIfstream fin;
        fin.openOrExit(inputFile);

        string ID;
        string line;
        while (fin >> ID) {
            mSampleIndexMap[ID] = mSampleIds.size();
            mSampleIds.push_back(ID);
            getline(fin, line);
            getline(fin, line); // ignore line for second haplotype
        }
        fin.close();

        int N = mSampleIds.size();
        int H = 2*N;
        if (mDebug > 0) {
            cout << "Read " << N << " IDs" << endl;
        }

        int Htop = 0;
        istringstream iss(line);
        while (iss >> line) {
            Htop++;
        }
        Htop /= 6;
        if (mDebug > 0) {
            cout << "Reading " << Htop << " IBS partners per line" << endl;
        }

        // read IBS lengths
        mIBSMatrix.resize(H, vector < pair <double, int> > (Htop));
        fin.openOrExit(inputFile);
        for (int h = 0; h < H; h++) {
            int hap12;
            fin >> ID >> hap12;
            assert(ID == mSampleIds[h/2]);
            for (int i = 0; i < Htop; i++) {
                double cMs[4];
                fin >> ID >> hap12 >> cMs[0] >> cMs[1] >> cMs[2] >> cMs[3];
                double lenLeft = weight0err * -cMs[1] + (1-weight0err) * -cMs[0];
                double lenRight = weight0err * cMs[2] + (1-weight0err) * cMs[3];
                double lenTot = lenLeft + lenRight;
                double lenMin = min(lenLeft, lenRight);
                mIBSMatrix[h][i].first = weightTotalLen * lenTot + (1-weightTotalLen) * lenMin;
                mIBSMatrix[h][i].second = 2*mSampleIndexMap[ID] + hap12 - 1;
                if (mIBSMatrix[h][i].second == -1) {
                    cout << "ERROR: " << h << " " << i << " " << ID << " " << hap12 << endl;
                    exit(1);
                }
            }
            sort(mIBSMatrix[h].begin(), mIBSMatrix[h].end(), greater < pair <double, int> > ());
        }
        fin.close();

        if (mVerbose > 0) {
            cout << timestamp() << " Read IBS matrix with " << mSampleIds.size() << " samples" << endl;
        }
    }

    vector< vector<double> > PhaseImpMissing::getDiploidCNPs(Variant* variant) {
        int H = mIBSMatrix.size();
        vector< vector<double> > dipCNPs(H/2);
        int numFound = 0;
        int numConfident = 0;
        vector< string > samples = variant->getSampleIds();
        vector< vector<float> > cnlVector = variant->getCNLs();
        for (uint i = 0; i < samples.size(); i++) {
            string sample = samples[i];
            vector<float>& cnls = cnlVector[i];
            vector<double> cnps = vector<double>(cnls.size());
            uint ncnls = cnls.size();
            double totalP = 0;
            double maxP = 0;
            for (uint j = 0; j < ncnls; j++) {
                cnps[j] = pow(10.0, cnls[j]);
                totalP += cnps[j];
            }
            for (uint j = 0; j < ncnls; j++) {
                cnps[j] = cnps[j] / totalP;
                if (cnps[j] > maxP) {
                    maxP = cnps[j];
                }
            }
            if (!mSampleIndexMap.count(sample)) {
                throw std::runtime_error(format("Sample not found in IBS map: %s", sample.c_str()));
            }
            int sampleIndex = mSampleIndexMap[sample];
            dipCNPs[sampleIndex] = cnps;
            numFound++;
            if (maxP >= 0.95) {
                numConfident++;
            }
        }
        if (mDebug > 0) {
            cout << "Read diploid CNLs for " << numFound << " individuals" << endl;
            cout << "(" << numConfident << " individuals have 95% confident dipCNs)" << endl;
        }
        return dipCNPs;
    }

    vector<string> PhaseImpMissing::getSampleStatusVector(const vector<string>& samples) {
        uint nSamples = samples.size();
        vector<string> PSTs(nSamples);
        for (uint i = 0; i < samples.size(); i++) {
            int sampleIndex = mSampleIndexMap[samples[i]];
            PSTs[i] = mSampleStatus[sampleIndex];
        }
        return PSTs;
    }

    static double computeCNF(const vector<double>& CNPs) {
        double cnf = 0;
        for (uint cn = 0; cn < CNPs.size(); cn++) {
            cnf += cn * CNPs[cn];
        }
        return cnf;
    }

    static double computeQual(vector<double>& CNPs) {
        // Temporarily overwrite CNPs to compute the second most likely value
        int cnMax1 = max_element(CNPs.begin(), CNPs.end()) - CNPs.begin();
        double cnProb1 = CNPs[cnMax1];
        CNPs[cnMax1] = 0;
        int cnMax2 = max_element(CNPs.begin(), CNPs.end()) - CNPs.begin();
        double cnProb2 = CNPs[cnMax2];
        CNPs[cnMax1] = cnProb1;
        double qual = 0.0;
        if (cnProb1 > cnProb2) {
            if (cnProb2 == 0) {
                qual = MAX_GENOTYPE_QUALITY;
            } else {
                qual = 10.0 * log10(cnProb1 / cnProb2);
                // cout << "#DBG: cnProb1=" << cnProb1 << " cnProb2=" << cnProb2 << " qual=" << qual << endl;
            }
            if (qual > MAX_GENOTYPE_QUALITY) {
                qual = MAX_GENOTYPE_QUALITY;
            }
        }
        return qual;
    }

    static string formatCNLs(const vector<double>& CNPs) {
        vector<double> CNLs = CNPs;
        uint outLength = 0;
        for (uint i = 0; i < CNLs.size(); i++) {
            CNLs[i] = log10(CNLs[i]);
            if (CNLs[i] < MIN_CN_LIKELIHOOD) {
                CNLs[i] = MIN_CN_LIKELIHOOD;
            } else {
                outLength = i;
            }
        }
        ostringstream buffer;
        for (uint i = 0; i < outLength; i++) {
            if (i > 0) {
                buffer << ",";
            }
            buffer << format("%1.1f", CNLs[i]);
        }
        return buffer.str();
    }

    static void formatSampleGenotypeAttrs(vector< vector<double> >& hapCNPs, string& PCN, string& PCNF, string& PCNQ, string& PCNL) {
        ostringstream buffer_PCN;
        ostringstream buffer_PCNF;
        ostringstream buffer_PCNQ;
        ostringstream buffer_PCNL;
        for (uint h = 0; h < hapCNPs.size(); h++) {
            int cnMax = max_element(hapCNPs[h].begin(), hapCNPs[h].end()) - hapCNPs[h].begin();
            double qual = computeQual(hapCNPs[h]);
            double cnf = computeCNF(hapCNPs[h]);
            if (h > 0) {
                buffer_PCN << "|";
                buffer_PCNF << "|";
                buffer_PCNQ << "|";
                buffer_PCNL << "|";
            }
            buffer_PCN << format("%d", cnMax);
            buffer_PCNF << format("%1.3f", cnf);
            buffer_PCNQ << format("%1.1f", qual);
            buffer_PCNL << formatCNLs(hapCNPs[h]);
        }
        PCN = buffer_PCN.str();
        PCNF = buffer_PCNF.str();
        PCNQ = buffer_PCNQ.str();
        PCNL = buffer_PCNL.str();
    }

    Variant* PhaseImpMissing::processVariant(Variant* variant) {
        vector< vector<double> > dipCNPs = getDiploidCNPs(variant);
        vector< vector<double> > hapCNPs = phaseImpCore(mIBSMatrix, dipCNPs, mIterations, mDebug);

        Variant* result = variant->reheader(mOutputHeader);
        vector<string> samples = result->getSampleIds();
        uint nSamples = samples.size();
        vector<string> PSTs = getSampleStatusVector(samples);
        vector<string> PCNs(nSamples);
        vector<string> PCNFs(nSamples);
        vector<string> PCNQs(nSamples);
        vector<string> PCNLs(nSamples);
        for (uint i = 0; i < samples.size(); i++) {
            int sampleIndex = mSampleIndexMap[samples[i]];
            vector< vector<double> > sampleHapCNPs(2);
            sampleHapCNPs[0] = hapCNPs[2*sampleIndex];
            sampleHapCNPs[1] = hapCNPs[2*sampleIndex + 1];
            formatSampleGenotypeAttrs(sampleHapCNPs, PCNs[i], PCNFs[i], PCNQs[i], PCNLs[i]);
        }

        result->updateFormatField("PST", PSTs);
        result->updateFormatField("PCN", PCNs);
        result->updateFormatField("PCNF", PCNFs);
        result->updateFormatField("PCNQ", PCNQs);
        result->updateFormatField("PCNL", PCNLs);
        return result;
    }
}
