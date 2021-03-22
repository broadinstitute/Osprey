
#ifndef OSPREY_PHASE_IMP_MISSING_HPP
#define OSPREY_PHASE_IMP_MISSING_HPP

#include <map>
#include <string>
#include <vector>

#include "Types.hpp"
#include "Variant.hpp"
#include "VCFReader.hpp"
#include "VCFWriter.hpp"
#include "vcf.h"

namespace Osprey {

    class PhaseImpMissing {
    private:
        int mDebug;
        int mVerbose;
        int mIterations;
        bcf_hdr_t* mOutputHeader;
        std::vector<std::string> mSampleIds;
        std::map<std::string, int> mSampleIndexMap;
        std::vector<std::string> mSampleStatus;
        std::vector< std::vector< std::pair<double, int> > > mIBSMatrix;

        std::vector< std::vector<double> > getDiploidCNPs(Variant* variant);
        std::vector<std::string> getSampleStatusVector(const std::vector<std::string>& samples);

    public:
        PhaseImpMissing();
        ~PhaseImpMissing();
        void setDebug(int value);
        void setVerbose(int value);
        void setIterations(int value);
        void readIBSMatrix(std::string inputFile);
        void updateVCFHeader(VCFReader& reader, VCFWriter& writer);
        Variant* processVariant(Variant* variant);
    };
}

#endif
