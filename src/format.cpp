
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "format.hpp"

std::string format(const std::string& format, ...)
{
    va_list args;
    va_start(args, format);
    size_t len = std::vsnprintf(NULL, 0, format.c_str(), args);
    va_end(args);
    std::vector<char> vec(len + 1);
    va_start(args, format);
    std::vsnprintf(&vec[0], len + 1, format.c_str(), args);
    va_end(args);
    return &vec[0];
}

std::string formatVector(const std::vector<int>& v)
{
    std::ostringstream str;
    str << "[";
    for (uint i = 0; i < v.size(); i++) {
        if (i > 0) {
            str << ",";
        }
        str << v[i];
    }
    str << "]";
    return(str.str());
}

std::string formatVector(const std::vector<uint>& v)
{
    std::ostringstream str;
    str << "[";
    for (uint i = 0; i < v.size(); i++) {
        if (i > 0) {
            str << ",";
        }
        str << v[i];
    }
    str << "]";
    return(str.str());
}
