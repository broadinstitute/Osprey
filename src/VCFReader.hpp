
#ifndef OSPREY_VCF_READER_HPP
#define OSPREY_VCF_READER_HPP

#include <string>
#include <vector>

#include "Types.hpp"
#include "Variant.hpp"
#include "vcf.h"

namespace Osprey {

    class VCFReader {
    private:
        std::string inputFile;
        htsFile* pBcfFile;
        bcf_hdr_t* pBcfHeader;

        void closeBcfFile();

    // internal
    public:
        bcf_hdr_t* getHeader();

    public:
        VCFReader(const std::string& inputFile);
        ~VCFReader();
        void close();
        std::vector<std::string> getSampleIds();
        Variant* nextVariant();
    };
}

#endif
