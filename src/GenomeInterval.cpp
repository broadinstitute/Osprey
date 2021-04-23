
#include <cstdlib>
#include <iostream>
#include <string>

#include "GenomeInterval.hpp"

using namespace std;
using namespace Osprey;

namespace Osprey {

    GenomeInterval::GenomeInterval() {
        start = 0;
        end = 0;
    }

    bool GenomeInterval::isNull() {
        return seqname.empty();
    }

    static bool parseCoordinate(const string& text, uint32_t& result) {
        char* eptr = NULL;
        long value = strtol(text.c_str(), &eptr, 10);
        if (*eptr != 0 || value < 0) {
            return false;
        }
        result = (uint32_t) value;
        return true;
    }

    bool GenomeInterval::parse(const string& text) {
        return parse(text.c_str());
    }

    bool GenomeInterval::parse(const char* text) {
        seqname.clear();
        start = 0;
        end = 0;
        if (text == NULL || text[0] == 0) {
            return false;
        }
        string input(text);
        size_t idx1 = input.find(':');
        if (idx1 == string::npos) {
            // Sequence name only is a valid interval
            seqname = input;
            return true;
        }
        size_t idx2 = input.find('-',idx1+1);
        if (idx2 == string::npos) {
            return false;
        }
        uint32_t startval;
        uint32_t endval;
        if (!parseCoordinate(input.substr(idx1+1,idx2-idx1-1), startval) ||
            !parseCoordinate(input.substr(idx2+1), endval)) {
            return false;
        }
        // Note we allow end == start-1 in order to represent zero-length intervals.
        if (endval + 1 < startval) {
            return false;
        }
        seqname = input.substr(0,idx1);
        start = startval;
        end = endval;
        return true;
    }
}

