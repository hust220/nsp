#include "Time.h"

namespace jian {

string Time::time() {
    time_t t;
    std::time(&t);
    string str = asctime(localtime(&t));
    return str.substr(0, str.size() - 1);
}

int Time::year() {
    return 0;
}

int Time::mon() {
    return 0;
}

int Time::date() {
    return 0;
}

int Time::hour() {
    int t = std::time(0);
    t = t % (3600 * 24);
    t = t / 3600 + 8;
    return t;
}

int Time::min() {
    int t = std::time(0);
    t = t % 3600;
    t = t / 60;
    return t;
}

int Time::sec() {
    int t = std::time(0);
    t = t % 60;
    return t;
}

} /// namespace jian
















