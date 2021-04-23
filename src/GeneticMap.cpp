
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "FileUtils.hpp"
#include "GeneticMap.hpp"
#include "format.hpp"

// Hack
#include "FileUtils.cpp"
#include "StringUtils.cpp"

using namespace std;
using namespace Osprey;

static const string EAGLE_GMAP_FILE_HEADER = "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)";
static const string SHAPEIT_GMAP_FILE_HEADER = "pos\tchr\tcM";

namespace Osprey {

    GeneticMap::GeneticMap() {
    }

    GeneticMap::GeneticMap(const string& path) {
        open(path);
    }

    static int split(const string& input, vector<string>& tokens, char delim = '\t') {
        string token;
        istringstream stream(input);
        while (getline(stream, token, delim)) {
            tokens.push_back(token);
        }
        return tokens.size();
    }

    void GeneticMap::open(const string& path) {
        mMap.clear();
        FileUtils::AutoGzIfstream fin;
        fin.openOrExit(path);
        string line;
        getline(fin, line);
        if (line == EAGLE_GMAP_FILE_HEADER) {
            string chr;
            uint32_t bp;
            double rate;
            double cm;
            while (fin >> chr >> bp >> rate >> cm) {
                if (!addEntry(chr, bp, cm)) {
                    throw std::runtime_error(format("Entries not sorted in genetic map file: %s", path.c_str()));
                }
            }
        } else if (line == SHAPEIT_GMAP_FILE_HEADER) {
            while (getline(fin, line)) {
                vector<string> tokens;
                if (split(line, tokens) != 3) {
                    throw std::runtime_error(format("Invalid line in genetic map file: %s", path.c_str()));
                }
                string chr;
                uint32_t bp;
                double cm;
                chr = tokens[1];
                sscanf(tokens[0].c_str(), "%u", &bp);
                sscanf(tokens[2].c_str(), "%lf", &cm);
                if (!addEntry(chr, bp, cm)) {
                    throw std::runtime_error(format("Entries not sorted in genetic map file: %s", path.c_str()));
                }
            }
        } else {
            throw std::runtime_error(format("Unrecognized file header in genetic map file: %s", path.c_str()));
        }
        fin.close();
    }

    vector< pair<uint32_t, double> >* GeneticMap::getSequenceMap(const string& chr) {
        typedef pair<uint32_t, double> PairType;
        typedef vector<PairType> VecType;
        typedef map<string, VecType> MapType;
        MapType::iterator iter = mMap.find(chr);
        // cout << "#DBG: getSequenceMap " << chr.c_str() << " iter=end " << (iter == mMap.end()) << endl;
        // cout << "#DBG: getSequenceMap " << chr.c_str() << " iter->key " << (iter->first.c_str()) << endl;
        if (iter == mMap.end() || iter->first != chr) {
            return NULL;
        }
        return &iter->second;
    }

    // This is forgiving about sequence names.
    vector< pair<uint32_t, double> >* GeneticMap::findSequenceMap(const string& chr) {
        vector< pair<uint32_t, double> >* result = getSequenceMap(chr);
        if (result == NULL) {
            // Try adding or removing a "chr" suffix to find a match.
            if (chr.length() > 3 && chr.substr(0,3) == "chr") {
                result = getSequenceMap(chr.substr(3));
            } else {
                result = getSequenceMap(format("chr%s", chr.c_str()));
            }
        }
        return result;
    }

    bool GeneticMap::addEntry(const string& chr, uint32_t pos, double cm) {
        typedef pair<uint32_t, double> PairType;
        typedef vector<PairType> VecType;
        typedef map<string, VecType> MapType;
        MapType::iterator iter = mMap.lower_bound(chr);
        if (iter == mMap.end() || iter->first != chr) {
            mMap.insert(MapType::value_type(chr, VecType()));
            iter = mMap.lower_bound(chr);
        }
        VecType& vec = iter->second;
        if (!vec.empty()) {
            PairType& lp = vec[vec.size()-1];
            if (lp.first >= pos || lp.second > cm) {
                return false;
            }
        }
        vec.push_back(PairType(pos, cm));
        return true;
    }

    struct PairComp {
        bool operator() (const pair<uint32_t, double>& pair, const uint32_t& key) {
            return (pair.first < key);
        }
        bool operator() (const uint32_t& key, const pair<uint32_t, double>& pair) {
            return (key < pair.first);
        }
    };

    double GeneticMap::interpolate(const string& chr, uint32_t pos) {
        typedef pair<uint32_t, double> PairType;
        typedef vector<PairType> VecType;
        VecType& vec = *findSequenceMap(chr);
        if (&vec == NULL) {
            throw std::runtime_error(format("Sequence not found in genetic map: %s", chr.c_str()));
        }
        if (vec.size() <= 1) {
            throw std::runtime_error(format("Cannot interpolate from sparse genetic map: %s:%d", chr.c_str(), pos));
        }
        VecType::const_iterator lbIter = lower_bound(vec.begin(), vec.end(), pos, PairComp());
        uint idx = lbIter - vec.begin();
        double result;
        if (idx == vec.size()) {
            double rate = (vec[idx-1].second - vec[idx-2].second) / (vec[idx-1].first - vec[idx-2].first);
            result = vec[idx-1].second + rate * (pos - vec[idx-1].first);
        } else if (vec[idx].first == pos) {
            result = vec[idx].second;
        } else if (idx == 0) {
            double rate = vec[idx].second / vec[idx].first;
            result = rate * pos;
        } else {
            double rate = (vec[idx].second - vec[idx-1].second) / (vec[idx].first - vec[idx-1].first);
            result = vec[idx-1].second + rate * (pos - vec[idx-1].first);
        }
        // cout << "#DBG: GeneticMap::interpolate " << chr << ":" << pos << " idx=" << idx << endl;
        // cout << "#DBG: vec size=" << vec.size() << endl;
        // cout << "#DBG: vec[" << idx << "] = " << vec[idx].first << "," << vec[idx].second << endl;
        // cout << "#DBG: result = " << result << endl;
        return result;
    }
}
