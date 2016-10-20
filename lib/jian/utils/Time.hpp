#pragma once

#include <string>
#include <ctime>

namespace jian {

class Time {
public:
    static std::string time() {
        time_t t;
        std::time(&t);
        std::string str = asctime(localtime(&t));
        return str.substr(0, str.size() - 1);
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

