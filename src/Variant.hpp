
#ifndef OSPREY_VARIANT_HPP
#define OSPREY_VARIANT_HPP

#include <string>
#include <vector>

#include "Types.hpp"
#include "vcf.h"

namespace Osprey {

    class Variant {
        friend class VCFWriter;

    private:
        bcf1_t* pRecord;
        bcf_hdr_t* pHeader;

        bcf1_t* getBcfRecord() { return pRecord; }

        std::vector< std::vector<float> > getFormatFloatVectors(const char* key);
        std::vector< std::vector<float> > unpackFmtFloatVectors(float* ptr, int nsamples, int nvalues);
        void unpack(int which);

    public:
        Variant(bcf1_t* record, bcf_hdr_t* header);
        ~Variant();

        std::string getChrom();
        uint32_t getPos();
        std::string getId();
        std::string getSiteText();
        std::vector<std::string> getSampleIds();
        std::vector<uint8_t> getEncodedGenotypes();
        std::vector< std::vector<float> > getCNLs();
        std::vector< std::vector<float> > getCNPs();

        void updateInfoField(const std::string& key, const std::string& value);
        void updateFormatField(const std::string& key, const std::vector<std::string>& values);
        Variant* reheader(bcf_hdr_t* header);
    };
}

#endif
