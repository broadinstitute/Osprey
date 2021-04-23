
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "omp.h"

#include "FileUtils.hpp"
#include "GeneticMap.hpp"
#include "GenomeInterval.hpp"
#include "OspreyIBSParams.hpp"
#include "Variant.hpp"
#include "VCFReader.hpp"
#include "format.hpp"
#include "timestamp.hpp"

using namespace std;
using namespace Osprey;

// Direct port of core algorithm for now.

static const int MAX_ERR = 2;
static const int MAX_MATCHES = 200;

struct IBSMatch {
    uint hap;
    uint hapKey;
    int32_t bpErrs[2][MAX_ERR];
    double cMerrs[2][MAX_ERR];
    double len;

    bool operator < (const IBSMatch& match2) const {
        if (len != match2.len) {
            return (len > match2.len);
        } else {
            // Randomized key for tie breaks
            return hapKey < match2.hapKey;
        }
    }
};

static vector<IBSMatch> findBestIBS(uint h1,
                                    uchar* genoPtrs[2],
                                    const vector<uint32_t> bps[2],
                                    const vector<double> cMs[2],
                                    const vector<uint> hapKeys,
                                    double weight0err,
                                    double weightTotalLen,
                                    const vector<bool>& refSampleMask) {

    const uint H = hapKeys.size();
    const uint flankCounts[2] = { (uint) bps[0].size(), (uint) bps[1].size() };

    set<IBSMatch> bestMatches;
    uint64_t num_checks = 0;
    for (uint h2 = 0; h2 < H; h2++) {
        if (h2 == h1) {
            continue;
        }
        if (!refSampleMask[h2/2]) {
            continue;
        }
        IBSMatch matchInfo;
        matchInfo.hap = h2;
        matchInfo.hapKey = hapKeys[h2];
        for (int k = 0; k < 2; k++) {
            if (flankCounts[k] == 0) {
                // at telomere, no flanking region
                for (int errs = 0; errs < MAX_ERR; errs++) {
                    matchInfo.bpErrs[k][errs] = 0;
                    matchInfo.cMerrs[k][errs] = 0;
                }
                continue;
            }
            int errs = 0;
            const uchar* genos_h1 = genoPtrs[k] + h1 * (uint64_t) flankCounts[k];
            const uchar* genos_h2 = genoPtrs[k] + h2 * (uint64_t) flankCounts[k];
            for (uint m = 0; m < flankCounts[k]; m++) {
                num_checks++;
                if (genos_h1[m] != genos_h2[m]) {
                    matchInfo.bpErrs[k][errs] = bps[k][m];
                    matchInfo.cMerrs[k][errs] = cMs[k][m];
                    errs++;
                    if (errs == MAX_ERR) {
                        break;
                    }
                }
            }
            while (errs < MAX_ERR) {
                matchInfo.bpErrs[k][errs] = bps[k].back();
                matchInfo.cMerrs[k][errs] = cMs[k].back();
                errs++;
            }
        }
        double lenLeft = weight0err*-matchInfo.cMerrs[0][0] + (1-weight0err)*-matchInfo.cMerrs[0][1];
        double lenRight = weight0err*matchInfo.cMerrs[1][0] + (1-weight0err)*matchInfo.cMerrs[1][1];
        double lenTot = lenLeft + lenRight;
        double lenMin = min(lenLeft, lenRight);
        matchInfo.len = weightTotalLen*lenTot + (1-weightTotalLen)*lenMin;
        bestMatches.insert(matchInfo);
        if (bestMatches.size() > MAX_MATCHES) {
            bestMatches.erase(--bestMatches.end());
        }
    }

    // cout << timestamp() << " Average checks " << h1 << ": " << (num_checks/(double)H) << endl;
    vector<IBSMatch> result(bestMatches.begin(), bestMatches.end());
    return result;
}

static void initGenotypeMatrix(uchar* gtMatrix, int N, int M, int m, const vector<uint8_t>& gts) {
    assert(gts.size() == (size_t) N);
    for (int n = 0; n < N; n++) {
        uint8_t gt = gts[n];
        if ((gt & 0x80) == 0) {
            throw std::runtime_error(format("Input genotypes are not phased."));
        }
        uchar a1 = (gt & 0x7);
        uchar a2 = ((gt >> 3) & 0x7);
        if (a1 > 1 || a2 > 1) {
            throw std::runtime_error(format("Input genotypes are not bi-allelic."));
        }
        gtMatrix[(2LL*n)*M + m] = a1;
        gtMatrix[(2LL*n+1)*M + m] = a2;
    }
}

static string formatGenotypes(const vector<uint8_t>& gts) {
    ostringstream oss;
    oss << "[";
    for (uint i = 0; i < gts.size(); i++) {
        if (i > 0) {
            oss << ",";
        }
        oss << format("%02X", gts[i]);
    }
    oss << "]";
    return oss.str();
}

int main(int argc, const char* argv[]) {

    const double weight0err = 0.5;
    const double weightTotalLen = 0.25;

    OspreyIBSParams params;
    GenomeInterval targetInterval;
    int status = params.processCommandLineArgs(argc, argv);
    if (status == 0) {
        exit(0);
    }
    if (status < 0) {
        // cerr << "Aborting due to error processing command line arguments" << endl;
        exit(-status);
    }

    if (params.threads > 64) {
        cerr << "Invalid number of threads: " << params.threads << endl;
        exit(1);
    }
    targetInterval.parse(params.ibsTargetInterval);
    if (targetInterval.isNull() || targetInterval.start == 0 || targetInterval.end == 0) {
        cerr << "Invalid target interval: " << params.ibsTargetInterval << endl;
        exit(1);
    }

    if (params.verbose > 0) {
        cout << timestamp() << " OspreyIBS version " << params.getVersionString() << endl;
        cout << timestamp() << " Input file: " << params.inputFile << endl;
        cout << timestamp() << " IBS matrix file: " << params.ibsMatrixFile << endl;
        cout << timestamp() << " Interval: " << params.ibsTargetInterval << endl;
        cout << timestamp() << " Genetic map file: " << params.geneticMapFile << endl;
        cout << timestamp() << " Reference samples: " << (params.refSampleList.empty() ? "all" : format("%d", params.refSampleList.size())) << endl;
        cout << timestamp() << " Threads: " << params.threads << endl;
    }

    const int seed = 1;
    srand(seed);
    omp_set_num_threads(params.threads);

    if (params.verbose > 0) {
        cout << timestamp() << " Reading genetic map file ..." << endl;
    }
    GeneticMap geneticMap(params.geneticMapFile);

    if (params.verbose > 0) {
        cout << timestamp() << " Reading input file (pass 1) ..." << endl;
    }

    int flankCounts[2] = {0, 0};
    vector<string> sampleIDs;
    vector<uint32_t> positions[2];

    VCFReader vcfReader(params.inputFile);
    while (true) {
        Variant* variant = vcfReader.nextVariant();
        if (variant == NULL) {
            break;
        }
        if (sampleIDs.empty()) {
            sampleIDs = variant->getSampleIds();
        }
        if (variant->getChrom() != targetInterval.seqname) {
            continue;
        }
        uint32_t pos = variant->getPos();
        if (pos < targetInterval.start) {
            flankCounts[0]++;
            positions[0].push_back(pos);
        } else if (pos > targetInterval.end) {
            flankCounts[1]++;
            positions[1].push_back(pos);
        } else {
            continue;
        }
    }
    vcfReader.close();

    // Reverse the left flank coordinates
    reverse(positions[0].begin(), positions[0].end());

    const uint N = sampleIDs.size();
    const uint H = 2 * N;

    uint nRefSamples = N;
    vector<bool> refSampleMask(N, params.refSampleList.empty());
    if (!params.refSampleList.empty()) {
        nRefSamples = 0;
        set<string> refSampleSet(params.refSampleList.begin(), params.refSampleList.end());
        for (uint i = 0; i < N; i++) {
            if (refSampleSet.find(sampleIDs[i]) != refSampleSet.end()) {
                refSampleMask[i] = true;
                nRefSamples++;
            }
        }
    }

    if (params.verbose > 0) {
        cout << timestamp() << " Found input data for " << sampleIDs.size() << " samples." << endl;
        cout << timestamp() << " Using " << nRefSamples << " samples as reference panel." << endl;
        cout << timestamp() << " Found " << flankCounts[0] << " variants"
             << " before " << targetInterval.seqname << ":" << targetInterval.start << "." << endl;
        cout << timestamp() << " Found " << flankCounts[1] << " variants"
             << " after " << targetInterval.seqname << ":" << targetInterval.end << "." << endl;
    }

    if (params.verbose > 0) {
        cout << timestamp() << " Reading input file (pass 2) ..." << endl;
    }

    vector<double> cmPositions[2];
    cmPositions[0].resize(flankCounts[0]);
    cmPositions[1].resize(flankCounts[1]);
    uint32_t midpos = targetInterval.start + (targetInterval.end - targetInterval.start) / 2;
    double cmMid = geneticMap.interpolate(targetInterval.seqname, midpos);
    int idxLeft = flankCounts[0] - 1;
    int idxRight = 0;

    uchar* genoPtrs[2];
    for (int k = 0; k < 2; k++) {
        size_t msize = H * (size_t) flankCounts[k];
        genoPtrs[k] = new uchar[msize];
        memset(genoPtrs[k], 0, msize);
    }

    VCFReader vcfReader2(params.inputFile);
    while (true) {
        Variant* variant = vcfReader2.nextVariant();
        if (variant == NULL) {
            break;
        }
        if (variant->getChrom() != targetInterval.seqname) {
            continue;
        }
        uint32_t pos = variant->getPos();
        if (pos >= targetInterval.start && pos <= targetInterval.end) {
            continue;
        }
        double cm = geneticMap.interpolate(variant->getChrom(), pos);
        vector<uint8_t> gts = variant->getEncodedGenotypes();
        if (params.debug > 0) {
            cout << "#DBG: variant " << variant->getSiteText() << " @" << variant->getChrom() << ":" << pos
                 << " " << format("%1.3f cM", cm) << " gts: " << formatGenotypes(gts) << endl;
        }
        if (pos < targetInterval.start) {
            cmPositions[0][idxLeft] = cm - cmMid;
            initGenotypeMatrix(genoPtrs[0], N, flankCounts[0], idxLeft, gts);
            idxLeft--;
        } else if (pos > targetInterval.end) {
            cmPositions[1][idxRight] = cm - cmMid;
            initGenotypeMatrix(genoPtrs[1], N, flankCounts[1], idxRight, gts);
            idxRight++;
        }
    }
    vcfReader2.close();

    if (params.verbose > 0) {
        cout << timestamp() << " Finished reading input file." << endl;
        cout << timestamp() << " Starting IBS computations ..." << endl;
    }

    vector<uint> hapKeys(H);
    for (uint h = 0; h < H; h++) {
        hapKeys[h] = h;
    }
    std::random_shuffle(hapKeys.begin(), hapKeys.end());

    vector< vector<IBSMatch> > matches(H);

#pragma omp parallel for
    for (uint h1 = 0; h1 < H; h1++) {
        matches[h1] = findBestIBS(h1, genoPtrs, positions, cmPositions, hapKeys, weight0err, weightTotalLen, refSampleMask);
    }

    if (params.verbose > 0) {
        cout << timestamp() << " Finished IBS computations." << endl;
        cout << timestamp() << " Wrinting output matrix ..." << endl;
    }

    FileUtils::AutoGzOfstream fout;
    fout.openOrExit(params.ibsMatrixFile);
    fout << std::setprecision(3) << std::fixed;
    for (uint h = 0; h < H; h++) {
        fout << sampleIDs[h/2] << "\t" << ((h%2)+1);
        for (uint m = 0; m < matches[h].size(); m++) {
            const IBSMatch& match = matches[h][m];
            fout << "\t" << sampleIDs[match.hap/2]
                 << "\t" << ((match.hap%2)+1)
                 << "\t" << match.cMerrs[0][1]
                 << "\t" << match.cMerrs[0][0]
                 << "\t" << match.cMerrs[1][0]
                 << "\t" << match.cMerrs[1][1];
        }
        fout << endl;
    }
    fout.close();

    if (params.verbose > 0) {
        cout << timestamp() << " Run complete." << endl;
    }
    return(0);
}

