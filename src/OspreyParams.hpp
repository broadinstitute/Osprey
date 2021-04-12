
#ifndef OSPREY_PARAMS_HPP
#define OSPREY_PARAMS_HPP

#include <string>
#include <vector>
#include "Types.hpp"

namespace Osprey {

    class OspreyParams {
    public:
        static const char* getVersionString();
    public:
        std::string inputFile;
        std::string outputFile;
        std::string ibsMatrixFile;
        std::vector<std::string> siteList;
        std::vector<uint> siteIndexRanges;
        std::string runBenchmark;
        std::vector<std::string> benchmarkSampleList;

        int debug;
        int verbose;
        int iterations;
        int threads;
        int benchmarkBatchSize;

        OspreyParams();
        int processCommandLineArgs(int argc, const char* argv[]);
    };
}

#endif
