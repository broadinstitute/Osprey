
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "OspreyIBSParams.hpp"
#include "format.hpp"
#include "version.hpp"

using namespace std;

namespace Osprey {

    static int DEFAULT_THREADS = 1;
    static int DEFAULT_DEBUG_LEVEL = 0;
    static int DEFAULT_VERBOSE_LEVEL = 2;

    OspreyIBSParams::OspreyIBSParams() {
        debug = DEFAULT_DEBUG_LEVEL;
        verbose = DEFAULT_VERBOSE_LEVEL;
        threads = DEFAULT_THREADS;
    }

    const char* OspreyIBSParams::getVersionString() {
        return OSPREY_VERSION;
    }

    static void printVersion(ostream& stream) {
        stream << "Osprey version " << OSPREY_VERSION << endl;
    }

    static void printUsage(ostream& stream) {
        printVersion(stream);
        stream << "Usage: ospreyIBS [options] inputFile" << endl;
        stream << "Options:" << endl;
        stream << " --inputFile [-vcf] file         Input VCF file." << endl;
        stream << " --ibsMatrixFile [-ibs] file.gz  Output IBS matrix file (required)." << endl;
        stream << " --ibsInterval [-L] interval     Target interval (required)." << endl;
        stream << " --geneticMapFile [-gmap] file   Genetic map file (required)." << endl;
        stream << " --refSampleFile [-rs] file      File listing reference samples (all input samples)." << endl;
        stream << " --threads [-t] N                Number of threads (" << DEFAULT_THREADS << ")." << endl;
        stream << " --help [-h]                     Print help synopsis." << endl;
        stream << " --version                       Print software version (" << OSPREY_VERSION << ")." << endl;
        stream << " --verbose N                     Verbosity level (" << DEFAULT_VERBOSE_LEVEL << ")." << endl;
        stream << " --debug N                       Debug output level (" << DEFAULT_DEBUG_LEVEL << ")." << endl;
    }

    static bool parseOption(int argc, const char** argv, int& argind, const char*& option, const char*& optarg) {
        if (argind >= argc) {
            return false;
        }
        optarg = NULL;
        const char* arg = argv[argind];
        if (arg != NULL && arg[0] == '-') {
            argind++;
            option = arg;
        } else {
            return false;
        }
        if (argind < argc) {
            arg = argv[argind];
            if (arg != NULL && arg[0] != '-') {
                argind++;
                optarg = arg;
            }
        }
        return true;
    }

    static bool matches(const char* str, const char* match) {
        return(str != NULL && match != NULL && !strcmp(str, match));
    }

    static bool requireArg(const char* option, const char* optarg) {
        if (optarg == NULL) {
            cerr << "Option " << option << " requires an argument" << endl;
            printUsage(cerr);
            return false;
        }
        return true;
    }

    static void argerror(const std::string& message) {
        cerr << "Error: " << message << endl;
        printUsage(cerr);
    }

    static string trimWhitespace(const std::string& str) {
        const std::string whitespace(" \t");
        const size_t strBegin = str.find_first_not_of(whitespace);
        if (strBegin == std::string::npos) {
            return "";
        }
        const size_t strEnd = str.find_last_not_of(whitespace);
        const size_t strRange = strEnd - strBegin + 1;
        return str.substr(strBegin, strRange);
    }

    static vector<string> parseStringListFromFile(const char* path) {
        string line;
        string token;
        vector<string> result;
        ifstream stream(path);
        while (getline(stream, line)) {
            istringstream st(line);
            if (getline(st, token, '\t')) {
                token = trimWhitespace(token);
                if (!token.empty()) {
                    result.push_back(token);
                }
            }
        }
        return result;
    }

    int OspreyIBSParams::processCommandLineArgs(int argc, const char* argv[]) {

        int argind = 1;
        const char* option = NULL;
        const char* optarg = NULL;

        while (parseOption(argc, argv, argind, option, optarg)) {
            if (matches(option, "-h") || matches(option, "--help")) {
                printUsage(cout);
                return 0;
            }
            if (matches(option, "--version")) {
                printVersion(cout);
                return 0;
            }
            if (matches(option, "-vcf") || matches(option, "--inputFile")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                if (optarg == NULL) {
                    argerror(format("Option %s requires an argument", option));
                    return -1;
                }
                inputFile = optarg;
            } else if (matches(option, "-ibs") || matches(option, "--ibsMatrixFile")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                ibsMatrixFile = optarg;
            } else if (matches(option, "--ibsPositionMatrixFile")) {
                // Hidden argument
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                ibsPositionMatrixFile = optarg;
            } else if (matches(option, "-L") || matches(option, "--ibsInterval")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                ibsTargetInterval = optarg;
            } else if (matches(option, "-gmap") || matches(option, "--geneticMapFile")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                geneticMapFile = optarg;
            } else if (matches(option, "-rs") || matches(option, "--refSampleFile")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                refSampleList = parseStringListFromFile(optarg);
            } else if (matches(option, "-t") || matches(option, "--threads")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                int value = atoi(optarg);
                if (value <= 0) {
                    argerror(format("Invalid value for %s: %s", option, optarg));
                    return -1;
                }
                threads = value;
            } else if (matches(option, "--verbose")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                int value = atoi(optarg);
                if (value < 0) {
                    value = 0;
                }
                verbose = value;
            } else if (matches(option, "--debug")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                int value = atoi(optarg);
                if (value < 0) {
                    value = 0;
                }
                debug = value;
            } else {
                argerror(format("Unrecognized option %s", option));
                return -1;
            }
        }

        int nPositional = argc - argind;
        if (inputFile.empty() && nPositional == 1) {
            inputFile = argv[argind];
        } else if (nPositional > 0) {
            if (!inputFile.empty()) {
                argerror("Input file must be specied as an option or positional argument but not both");
            } else {
                argerror("Too many arguments supplied");
            }
            return -1;
        }

        if (inputFile.empty()) {
            argerror("Input file argument is required");
            return -1;
        }
        if (ibsMatrixFile.empty()) {
            argerror("Option --ibsMatrixFile [-ibs] is required");
            return -1;
        }
        if (ibsTargetInterval.empty()) {
            argerror("Option --ibsInterval [-L] is required");
            return -1;
        }
        if (geneticMapFile.empty()) {
            argerror("Option --geneticMapFile [-gmap] is required");
            return -1;
        }

        // cout << "Input file: " << inputFile << endl;
        // cout << "IBS matrix file: " << ibsMatrixFile << endl;
        // cout << "Interval: " << ibsTargetInterval << endl;
        // cout << "Genetic map file: " << geneticMapFile << endl;
        // cout << "Reference sample file: " << refSampleFile << endl;
        // cout << "Threads: " << threads << endl;

        return 1;
    }
}

  
