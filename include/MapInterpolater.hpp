#ifndef MAPINTERPOLATER_HPP
#define MAPINTERPOLATER_HPP

#include <string>
#include <map>
#include <utility>

namespace Genetics {

  class MapInterpolater {
    std::map < std::pair <int, int>, std::pair <double, double> > chrBpToRateGen;
    static const std::string MAP_FILE_HEADER;
  public:
    // input file format: chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
    // (Oxford map format preceded by chr column)
    MapInterpolater(const std::string &geneticMapFile);
    // returns interpolated genetic position in Morgans
    double interp(int chr, int bp) const;
  };

}
#endif
