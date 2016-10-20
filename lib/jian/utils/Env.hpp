#pragma once

#include <string>
#include <iostream>

namespace jian {

class Env {
public:
    std::string _lib;

    static Env &instance() {
        static Env env;
        return env;
    }

    static std::string lib(const std::string &s = "") {
        if (s.empty()) {
            if (instance()._lib.empty()) {
                char *path = std::getenv("NSP");
                if (path == NULL) {
                    std::cout << "Please tell me the path of the NSP library: ";
                    std::cin >> instance()._lib;
                } else instance()._lib = path;
            }
        } else instance()._lib = s;
        return instance()._lib;
    }
};

}
