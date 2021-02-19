#ifndef NUMERICUTILS_HPP
#define NUMERICUTILS_HPP

#include <cstdlib>
#include <vector>
#include <utility>

#include "Types.hpp"

namespace NumericUtils {

  inline double sq(double x) { return x*x; }
  double sum(const double x[], uint64 N);
  double mean(const std::vector <double> &x);

  // takes into account that some 0 values may indicate missing/ignored: divide out by Nused, not N
  double mean(const double x[], uint64 N, uint64 Nused);

  // regress y on x, assuming both have been 0-centered (so 0-filled missing values ok)
  double regCoeff(const double y[], const double x[], uint64 N);

  double dot(const double x[], const double y[], uint64 N);
  double norm2(const double x[], uint64 N);
  void normalize(double x[], uint64 N);

  void logSumExp(float &x, float y);

  std::pair <double, double> meanStdDev(const double x[], uint64 N);
  std::pair <double, double> meanStdErr(const double x[], uint64 N);
  std::pair <double, double> meanStdDev(const std::vector <double> &x);
  std::pair <double, double> meanStdErr(const std::vector <double> &x);
}

#endif
