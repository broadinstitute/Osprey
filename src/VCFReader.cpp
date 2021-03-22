
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "vcf.h"

#include "VCFReader.hpp"
#include "format.hpp"

using namespace std;
using namespace Osprey;

namespace Osprey {

    VCFReader::VCFReader(const std::string& inputFile) {
        this->inputFile = inputFile;
        pBcfFile = bcf_open(inputFile.c_str(), "r");
        if (pBcfFile == NULL) {
            throw std::runtime_error(format("Unable to read input file: %s", inputFile.c_str()));
        }
        pBcfHeader = bcf_hdr_read(pBcfFile);
        if (pBcfHeader == NULL) {
            throw std::runtime_error(format("Unable to read input file header: %s", inputFile.c_str()));
        }
    }
    
    VCFReader::~VCFReader() {
        close();
    }

    void VCFReader::close() {
        if (pBcfHeader != NULL) {
            bcf_hdr_destroy(pBcfHeader);
            pBcfHeader = NULL;
        }
        if (pBcfFile != NULL) {
            bcf_close(pBcfFile);
            pBcfFile = NULL;
        }
    }

    bcf_hdr_t* VCFReader::getHeader() {
        return pBcfHeader;
    }

    vector<string> VCFReader::getSampleIds() {
        int nSamples = bcf_hdr_nsamples(pBcfHeader);
        vector<string> result(nSamples);
        for (int i = 0; i < nSamples; i++) {
            result[i] = pBcfHeader->samples[i];
        }
        return result;
    }

    Variant* VCFReader::nextVariant() {
        if (pBcfFile == NULL || pBcfHeader == NULL) {
            return NULL;
        }
        bcf1_t* pRecord = bcf_init();
        int status = bcf_read(pBcfFile, pBcfHeader, pRecord);
        if (status == -1) {
            // end of file
            close();
            return NULL;
        }
        if (status < 0) {
            throw std::runtime_error(format("Error reading input file %s: bcf_read status code %d", inputFile, status));
        }
        if (pRecord->errcode != 0) {
            throw std::runtime_error(format("Error reading input file %s: bcf_read error code %d", inputFile, pRecord->errcode));
        }
        return new Variant(pRecord, pBcfHeader);
    }

}
