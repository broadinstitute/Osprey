
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "vcf.h"

#include "Variant.hpp"
#include "format.hpp"

using namespace std;

namespace Osprey {

    // BCF2 MISSING marker
    /** Currently unused
    static bool isMissing(float value) {
        uint8_t* p = (uint8_t*) &value;
        uint32_t u = *(uint32_t*) p;
        return (u == 0x7F800001);
    }
    **/

    // BCF2 END_OF_VECTOR marker
    static bool isEOV(float value) {
        uint8_t* p = (uint8_t*) &value;
        uint32_t u = *(uint32_t*) p;
        return (u == 0x7F800002);
    }

    Variant::Variant(bcf1_t* record, bcf_hdr_t* header) {
        pRecord = record;
        pHeader = header;
    }

    Variant::~Variant() {
        bcf_destroy(pRecord);
        pRecord = NULL;
        pHeader = NULL;
    }

    std::string Variant::getId() {
        unpack(BCF_UN_STR);
        std::string id = pRecord->d.id;
        return id;
    }

    std::string Variant::getSiteText() {
        // Temporary
        return getId();
    }

    vector<string> Variant::getSampleIds() {
        int nSamples = bcf_hdr_nsamples(pHeader);
        vector<string> result(nSamples);
        for (int i = 0; i < nSamples; i++) {
            result[i] = pHeader->samples[i];
        }
        return result;
    }

    vector< vector<float> > Variant::getCNLs() {
        unpack(BCF_UN_FMT);
        bcf_fmt_t* fmt = bcf_get_fmt(pHeader, pRecord, "CNL");
        if (fmt == NULL) {
            return vector< vector<float> >();
        }
        if (fmt->type != BCF_BT_FLOAT) {
            throw std::runtime_error(format("Invalid type for CNL field: %d", fmt->type));
        }
        int nsamples = bcf_hdr_nsamples(pHeader);
        return unpackFmtFloatVectors((float*) fmt->p, nsamples, fmt->n);
    }

    vector< vector<float> > Variant::unpackFmtFloatVectors(float* ptr, int nsamples, int nvalues) {
        vector< vector<float> > result(nsamples, vector<float>());
        int index = 0;
        for (int i = 0; i < nsamples; i++) {
            for (int j = 0; j < nvalues; j++) {
                float value = ptr[index++];
                if (!isEOV(value)) {
                    result[i].push_back(value);
                }
            }
        }
        return result;
    }

    void Variant::updateFormatField(const string& key, const vector<string>& values) {
        const int count = (int) values.size();
        const char* buffer[count];
        for (int i = 0; i < count; i++) {
            buffer[i] = values[i].c_str();
        }
        int status = bcf_update_format_string(pHeader, pRecord, key.c_str(), buffer, count);
        if (status != 0) {
            throw std::runtime_error(format("Error from bcf_update_format: status code %d", status));
        }
    }

    static void reheaderConvertFormatValues(bcf_hdr_t* header, bcf_fmt_t* fmt, const char* tag, int* imap, void** values, int* nvalues) {
        const int nsamples = bcf_hdr_nsamples(header);
        const int n = fmt->n;
        const void* oldvalues = *values;
        void* newvalues = NULL;
        switch (fmt->type) {
            case BCF_BT_INT8:
            case BCF_BT_INT16:
            case BCF_BT_INT32:
                newvalues = malloc(sizeof(uint32_t) * n * nsamples);
                for (int i = 0; i < nsamples; i++) {
                    int idx = imap[i];
                    uint32_t missing_value = bcf_int32_missing;
                    uint32_t vector_end_value = bcf_int32_vector_end;
                    if (!strcmp(tag, "GT")) {
                        missing_value = bcf_gt_missing;
                    }
                    if (idx < 0) {
                        for (int j = 0; j < n; j++) {
                            ((uint32_t*)newvalues)[i*n+j] = (j == 0) ? missing_value : vector_end_value;
                        }
                    } else {
                        for (int j = 0; j < n; j++) {
                            ((uint32_t*)newvalues)[i*n+j] = ((uint32_t*)oldvalues)[idx*n+j];
                        }
                    }
                }
                break;
            case BCF_BT_FLOAT:
                newvalues = malloc(sizeof(float) * n * nsamples);
                for (int i = 0; i < nsamples; i++) {
                    int idx = imap[i];
                    uint32_t missing_value = bcf_float_missing;
                    uint32_t vector_end_value = bcf_float_vector_end;
                    if (idx < 0) {
                        for (int j = 0; j < n; j++) {
                            bcf_float_set(&((float*)newvalues)[i*n+j], (j == 0) ? missing_value : vector_end_value);
                        }
                    } else {
                        for (int j = 0; j < n; j++) {
                            bcf_float_set(&((float*)newvalues)[i*n+j], ((uint32_t*)oldvalues)[idx*n+j]);
                        }
                    }
                }
                break;
            case BCF_BT_CHAR:
                newvalues = malloc(sizeof(uint8_t) * n * nsamples);
                for (int i = 0; i < nsamples; i++) {
                    int idx = imap[i];
                    uint8_t missing_value = bcf_str_missing;
                    uint8_t vector_end_value = bcf_str_vector_end;
                    if (idx < 0) {
                        for (int j = 0; j < n; j++) {
                            ((uint8_t*)newvalues)[i*n+j] = (j == 0) ? missing_value : vector_end_value;
                        }
                    } else {
                        for (int j = 0; j < n; j++) {
                            ((uint8_t*)newvalues)[i*n+j] = ((uint8_t*)oldvalues)[idx*n+j];
                        }
                    }
                }
                break;
            case BCF_BT_INT64:
                // Unsupported: It is not clear how to set such values through the bcf_update_format API.
            default:
                throw std::runtime_error(format("Unsupported format type: %d", fmt->type));
        }
        free(*values);
        *values = newvalues;
        *nvalues = (n * nsamples);
    }

    Variant* Variant::reheader(bcf_hdr_t* newHeader) {

        bcf1_t* oldRecord = pRecord;
        bcf1_t* newRecord = bcf_init();
        bcf_hdr_t* oldHeader = pHeader;

        unpack(BCF_UN_ALL);
        const char* site = oldRecord->d.id;

        // Note: The code below assumes in some cases that the header is compatible (e.g. the FILTER IDs).

        newRecord->pos = oldRecord->pos;
        newRecord->rlen = oldRecord->rlen;
        newRecord->rid = oldRecord->rid;
        newRecord->qual = oldRecord->qual;

        int status = bcf_update_id(newHeader, newRecord, site);
        if (status != 0) {
            throw std::runtime_error(format("Error from bcf_update_id: site %s, status code %d", site, status));
        }
        status = bcf_update_filter(newHeader, newRecord, oldRecord->d.flt, oldRecord->d.n_flt);
        if (status != 0) {
            throw std::runtime_error(format("Error from bcf_update_filter: site %s, status code %d", site, status));
        }
        status = bcf_update_alleles(newHeader, newRecord, (const char**) oldRecord->d.allele, oldRecord->n_allele);
        if (status != 0) {
            throw std::runtime_error(format("Error from bcf_update_alleles: site %s, status code %d", site, status));
        }

        int nInfo = oldRecord->n_info;
        for (int i = 0; i < nInfo; i++) {
            bcf_info_t* inf = &oldRecord->d.info[i];
            const char* tag = bcf_hdr_int2id(oldHeader, BCF_DT_ID, inf->key);
            int htype = bcf_hdr_id2type(oldHeader, BCF_HL_INFO, inf->key);
            void* values = NULL;
            int nvalues = 0;
            // cout << format("#DBG: site %s, tag %s INFO", site, tag) << " htype " << htype << endl;
            status = bcf_get_info_values(oldHeader, oldRecord, tag, &values, &nvalues, htype);
            if (status < 0) {
                throw std::runtime_error(format("Error from bcf_get_info_values: site %s, tag %s, status code %d", site, tag, status));
            }
            status = bcf_update_info(newHeader, newRecord, tag, values, nvalues, htype);
            if (status < 0) {
                throw std::runtime_error(format("Error from bcf_update_info: site %s, tag %s, status code %d", site, tag, status));
            }
        }

        // For format values, for compute sample map.
        // sampleMap[i] is index of source sample (or -1 if sample is no source sample)
        bool samplesDiffer = 0;
        int nFormat = oldRecord->n_fmt;
        int nNewSamples = bcf_hdr_nsamples(newHeader);
        vector<int> sampleMap(nNewSamples, -1);
        if (nFormat > 0) {
            for (int i = 0; i < nNewSamples; i++) {
                sampleMap[i] = bcf_hdr_id2int(oldHeader, BCF_DT_SAMPLE, newHeader->samples[i]);
                if (sampleMap[i] != i) {
                    samplesDiffer = 1;
                }
            }
        }
        for (int i = 0; i < nFormat; i++) {
            bcf_fmt_t* fmt = &oldRecord->d.fmt[i];
            const char* tag = bcf_hdr_int2id(oldHeader, BCF_DT_ID, fmt->id);
            int newId = bcf_hdr_id2int(newHeader, BCF_DT_ID, tag);
            if (newId < 0) {
                // Field is dropped in new header
                continue;
            }
            int htype = bcf_hdr_id2type(oldHeader, BCF_HL_FMT, fmt->id);
            void* values = NULL;
            int nvalues = 0;
            int status = 0;
            // cout << format("#DBG: site %s, tag %s FMT", site, tag) << " fmt->type " << fmt->type << " htype " << htype << endl;
            if (!strcmp(tag, "GT")) {
                status = bcf_get_genotypes(oldHeader, oldRecord, &values, &nvalues);
            } else {
                status = bcf_get_format_values(oldHeader, oldRecord, tag, &values, &nvalues, htype);
            }
            if (status < 0) {
                throw std::runtime_error(format("Error from bcf_get_format_values: site %s, tag %s, status code %d", site, tag, status));
            }
            if (samplesDiffer) {
                // cout << "#DBG: before convert valp=" << format("%p", values) << " nvalues=" << nvalues << endl;
                reheaderConvertFormatValues(newHeader, fmt, tag, sampleMap.data(), &values, &nvalues);
                // cout << "#DBG: after convert valp=" << format("%p", values) << " nvalues=" << nvalues << endl;
            }
            if (!strcmp(tag,"GT")) {
                status = bcf_update_genotypes(newHeader, newRecord, values, nvalues);
            } else {
                status = bcf_update_format(newHeader, newRecord, tag, values, nvalues, htype);
            }
            if (status < 0) {
                throw std::runtime_error(format("Error from bcf_update_format site %s, tag %s, status code %d", site, tag, status));
            }
            free(values);
        }

        return new Variant(newRecord, newHeader);
    }

    void Variant::unpack(int which) {
        int status = bcf_unpack(pRecord, which);
        if (status != 0) {
            throw std::runtime_error(format("Error reading input file: bcf_unpack error code %d", status));
        }
    }
}
