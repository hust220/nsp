#ifndef JIAN_TIME_H
#define JIAN_TIME_H

#include "std.h"

namespace jian {

class Time {
public:
    static std::string time() {
        time_t t;
        std::time(&t);
        string str = asctime(localtime(&t));
        return str.substr(0, str.size() - 1);
    }

    static int year() {
        return 0;
    }

    static int mon() {
        return 0;
    }

    static int date() {
        return 0;
    }

    static int hour() {
        int t = std::time(0);
        t = t % (3600 * 24);
        t = t / 3600 + 8;
        return t;
    }

    static int min() {
        int t = std::time(0);
        t = t % 3600;
        t = t / 60;
        return t;
    }

    static int sec() {
        int t = std::time(0);
        t = t % 60;
        return t;
    }

};

}

#endif

