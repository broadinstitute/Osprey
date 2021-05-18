
#ifndef OSPREY_IBS_PARAMS_HPP
#define OSPREY_IBS_PARAMS_HPP

#include <string>
#include <vector>
#include "Types.hpp"

namespace Osprey {

    class OspreyIBSParams {
    public:
        static const char* getVersionString();
    public:
        std::string inputFile;
        std::string ibsMatrixFile;
        std::string ibsPositionMatrixFile;
        std::string ibsTargetInterval;
        std::string geneticMapFile;

        std::vector<std::string> refSampleList;

        int debug;
        int verbose;
        int threads;

        OspreyIBSParams();
        int processCommandLineArgs(int argc, const char* argv[]);
    };
}

#endif
