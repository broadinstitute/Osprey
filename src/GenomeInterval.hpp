
#ifndef OSPREY_GENOME_INTERVAL_HPP
#define OSPREY_GENOME_INTERVAL_HPP

#include <string>
#include "Types.hpp"

namespace Osprey {

    class GenomeInterval {
    public:
        std::string seqname;
        uint32_t start;
        uint32_t end;
    public:
        GenomeInterval();
        bool parse(const std::string& text);
        bool parse(const char* text);
        bool isNull();
    };
}

#endif
