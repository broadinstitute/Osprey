
#ifndef OSPREY_VCF_WRITER_HPP
#define OSPREY_VCF_WRITER_HPP

#include <string>
#include <vector>

#include "Types.hpp"
#include "Variant.hpp"
#include "vcf.h"

namespace Osprey {

    class VCFWriter {
    private:
        std::string outputFile;
        htsFile* pBcfFile;
        bcf_hdr_t* pBcfHeader;

    public:
        VCFWriter(const std::string& outputFile);
        ~VCFWriter();
        void close();
        bcf_hdr_t* getHeader();
        void writeHeader(bcf_hdr_t* header);
        void writeVariant(Variant* variant);
    };
}

#endif
