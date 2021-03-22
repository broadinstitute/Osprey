
#include <cstdlib>
#include <ctime>

#include "timestamp.hpp"

using namespace std;

string timestamp() {
    struct tm tm;
    char buffer[100];
    const time_t now = time(NULL);
    localtime_r(&now, &tm);
    strftime(buffer, sizeof(buffer), "%c", &tm);
    return string(buffer);
}
