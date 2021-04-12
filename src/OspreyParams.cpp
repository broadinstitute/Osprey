
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <fstream>
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
        stream << " --site ID                       Specify site to process or space separated list (all)." << endl;
        stream << " --siteIndex N@B | N[-M],...     Specify sites to process by variant index ranges or blocks (all)." << endl;
        stream << " --benchmark type                Run internal benchmarks (cross)." << endl;
        stream << " --benchmarkBatchSize N          Benchmark match size (1)." << endl;
        stream << " --benchmarkSamples file         List of specific samples to use for benchmarking." << endl;
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

    static vector<string> parseStringList(const char* arg) {
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
                    // Ideally we should strip leading/trailing whitespace
                    result.push_back(token);
                }
            }
        }
        return result;
    }

    static bool parseSiteIndex(const string& token, uint& value) {
        const char* p = token.c_str();
        char* end = NULL;
        long val = strtol(p, &end, 10);
        if (token.length() != (uint) (end - p)) {
            return false;
        }
        if (val < 1 || val > INT_MAX) {
            return false;
        }
        value = (uint) val;
        return true;
    }

    static bool parseOneIndex(const string& token, vector<uint>& values) {
        uint val;
        if (!parseSiteIndex(token, val)) {
            return false;
        }
        values.push_back(val);
        return true;
    }

    static bool parseTwoIndexes(const string& token, vector<uint>& values, const char sep) {
        size_t idx = token.find(sep);
        if (idx == string::npos || idx == 0 || idx >= token.size()-1) {
            return false;
        }
        if (token.find(sep, idx+1) != string::npos) {
            return false;
        }
        uint val1;
        uint val2;
        if (!parseSiteIndex(token.substr(0,idx), val1) ||
            !parseSiteIndex(token.substr(idx+1), val2)) {
            return false;
        }
        values.push_back(val1);
        values.push_back(val2);
        return true;
    }

    static bool mergeSiteRanges(vector<uint>& ranges, const vector<uint>& values) {
        assert(values.size() == 2);
        assert(values[0] > 0);
        assert(values[1] >= values[0]);
        if (ranges.empty()) {
            ranges = values;
        } else if (ranges[ranges.size()-1] == values[0]-1) {
            ranges[ranges.size()-1] = values[1];
        } else if (ranges[ranges.size()-1] < values[0]) {
            ranges.insert(ranges.end(), values.begin(), values.end());
        } else {
            return false;
        }
        return true;
    }

    static bool parseSiteIndexSpec(const char* arg, vector<uint>& ranges) {
        string str(arg);
        istringstream stream(str);
        string token;
        vector<uint> values;
        while (getline(stream, token, ',')) {
            values.clear();
            if (parseTwoIndexes(token, values, '@')) {
                uint n = values[0];
                uint b = values[1];
                values[0] = (n*b - b + 1);
                values[1] = (n*b);
                if (!mergeSiteRanges(ranges, values)) {
                    argerror(format("Invalid or conflicting site index specification: %s", arg));
                    return false;
                }
            } else if (parseTwoIndexes(token, values, '-')) {
                if (!mergeSiteRanges(ranges, values)) {
                    argerror(format("Invalid or conflicting site index specification: %s", arg));
                    return false;
                }
            } else if (parseOneIndex(token, values)) {
                values.push_back(values[0]);
                if (!mergeSiteRanges(ranges, values)) {
                    argerror(format("Invalid or conflicting site index specification: %s", arg));
                    return false;
                }
            } else {
                argerror(format("Invalid site index specification: %s", arg));
                return false;
            }
        }
        return true;
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
                vector<string> argvals = parseStringList(optarg);
                siteList.insert(siteList.end(), argvals.begin(), argvals.end());
            } else if (matches(option, "--siteIndex")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                if (!parseSiteIndexSpec(optarg, siteIndexRanges)) {
                    return -1;
                }
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
            } else if (matches(option, "--benchmark")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                if (matches(optarg, "cross")) {
                    runBenchmark = optarg;
                } else {
                    argerror(format("Unrecognized benchmark type: %s", optarg));
                    return -1;
                }
            } else if (matches(option, "--benchmarkBatchSize")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                int value = atoi(optarg);
                if (value <= 0) {
                    argerror(format("Invalid value for %s: %s", option, optarg));
                    return -1;
                }
                benchmarkBatchSize = value;
            } else if (matches(option, "--benchmarkSamples")) {
                if (!requireArg(option, optarg)) {
                    return -1;
                }
                benchmarkSampleList = parseStringListFromFile(optarg);
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

  
