
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "OspreyParams.hpp"
#include "format.hpp"
#include "version.hpp"

using namespace std;

namespace Osprey {

    static int DEFAULT_THREADS = 1;
    static int DEFAULT_ITERATIONS = 250;
    static int DEFAULT_DEBUG_LEVEL = 0;
    static int DEFAULT_VERBOSE_LEVEL = 2;

    OspreyParams::OspreyParams() {
        debug = DEFAULT_DEBUG_LEVEL;
        verbose = DEFAULT_VERBOSE_LEVEL;
        threads = DEFAULT_THREADS;
        iterations = DEFAULT_ITERATIONS;
    }

    const char* OspreyParams::getVersionString() {
        return OSPREY_VERSION;
    }

    static void printUsage(ostream& stream) {
        stream << "Usage: osprey [options] inputFile" << endl;
        stream << "Options:" << endl;
        stream << " --inputFile [-vcf] file         Input VCF file." << endl;
        stream << " --outputFile [-o] file          Output VCF file (required)." << endl;
        stream << " --ibsMatrixFile [-ibs] file     IBS matrix file (required)." << endl;
        stream << " --threads [-t] N                Number of threads (" << DEFAULT_THREADS << ")." << endl;
        stream << " --iterations [-iter] N          Number of phasing iterations (" << DEFAULT_ITERATIONS << ")." << endl;
        stream << " --site ID                       Specific site to process or space separated list." << endl;
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

    static vector<string> parse_string_list(const char* arg) {
        string str(arg);
        istringstream stream(str);
        string token;
        vector<string> result;
        while (getline(stream, token, ' ')) {
            if (!token.empty()) {
                result.push_back(token);
            }
        }
        return result;
    }

    int OspreyParams::processCommandLineArgs(int argc, const char* argv[]) {

        int argind = 1;
        const char* option = NULL;
        const char* optarg = NULL;

        while (parseOption(argc, argv, argind, option, optarg)) {
            if (matches(option, "-h") || matches(option, "--help")) {
                printUsage(cout);
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
            } else if (matches(option, "-o") || matches(option, "--outputFile")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                outputFile = optarg;
            } else if (matches(option, "-ibs") || matches(option, "--ibsMatrixFile")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                ibsMatrixFile = optarg;
            } else if (matches(option, "--site")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                vector<string> argvals = parse_string_list(optarg);
                siteList.insert(siteList.end(), argvals.begin(), argvals.end());
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
            } else if (matches(option, "-iter") || matches(option, "--iterations")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                int value = atoi(optarg);
                if (value <= 0) {
                    argerror(format("Invalid value for %s: %s", option, optarg));
                    return -1;
                }
                iterations = value;
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
        if (outputFile.empty()) {
            argerror("Option --outputFile [-o] is required");
            return -1;
        }
        if (ibsMatrixFile.empty()) {
            argerror("Option --ibsMatrixFile [-ibs] is required");
            return -1;
        }

        // cout << "Input file: " << inputFile << endl;
        // cout << "Output file: " << outputFile << endl;
        // cout << "IBS matrix file: " << ibsMatrixFile << endl;
        // cout << "Iterations: " << iterations << endl;
        // cout << "Threads: " << threads << endl;

        return 1;
    }
}

  
