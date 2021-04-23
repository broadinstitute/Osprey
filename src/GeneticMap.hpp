
#ifndef OSPREY_GENETIC_MAP_HPP
#define OSPREY_GENETIC_MAP_HPP

#include <map>
#include <string>
#include <utility>
#include <vector>
#include "Types.hpp"

namespace Osprey {

    class GeneticMap {
    private:
        std::map<std::string, std::vector< std::pair<uint32_t, double> > > mMap;

    private:
        std::vector< std::pair<uint32_t, double> >* getSequenceMap(const std::string& chr);
        std::vector< std::pair<uint32_t, double> >* findSequenceMap(const std::string& chr);
        bool addEntry(const std::string& seqname, uint32_t pos, double cm);

    public:
        GeneticMap();
        GeneticMap(const std::string& path);
        void open(const std::string& path);

        // Interpolated value from map in centimorgans.
        // Throws an exception if sequence name is not in the map (or if map is too sparse to interpolate).
        double interpolate(const std::string& chr, uint32_t pos);
    };
}

#endif
