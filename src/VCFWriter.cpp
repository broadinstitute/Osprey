
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "vcf.h"

#include "VCFWriter.hpp"
#include "format.hpp"

using namespace std;
using namespace Osprey;

namespace Osprey {

    VCFWriter::VCFWriter(const std::string& outputFile) {
        this->outputFile = outputFile;
        pBcfHeader = NULL;
        pBcfFile = bcf_open(outputFile.c_str(), "wz");
        if (pBcfFile == NULL) {
            throw std::runtime_error(format("Unable to open output file %s: %s", outputFile.c_str(), strerror(errno)));
        }
    }

    VCFWriter::~VCFWriter() {
        close();
    }

    void VCFWriter::close() {
        if (pBcfFile != NULL) {
            bcf_close(pBcfFile);
            pBcfFile = NULL;
        }
        if (pBcfHeader != NULL) {
            bcf_hdr_destroy(pBcfHeader);
            pBcfHeader = NULL;
        }
    }

    bcf_hdr_t* VCFWriter::getHeader() {
        return pBcfHeader;
    }

    void VCFWriter::writeHeader(bcf_hdr_t* header) {
        if (pBcfHeader != NULL) {
            throw std::runtime_error(format("Attempt to write output file header more than once: %s", outputFile.c_str()));
        }
        int status = bcf_hdr_write(pBcfFile, header);
        if (status != 0) {
            throw std::runtime_error(format("Error writing output file header for %s: status %d", outputFile.c_str(), status));
        }
        pBcfHeader = header;
    }

    void VCFWriter::writeVariant(Variant* variant) {
        if (pBcfFile == NULL || pBcfHeader == NULL) {
            throw std::runtime_error(format("Error writing output file %s: File is closed or no header written", outputFile.c_str()));
        }
        int status = bcf_write(pBcfFile, pBcfHeader, variant->getBcfRecord());
        if (status != 0) {
            throw std::runtime_error(format("Error writing output file record %s: status %d", outputFile.c_str(), status));
        }
    }
}
