#ifndef JIAN_UTIL_LIB
#define JIAN_UTIL_LIB

#include <string>
#include "../lib.hpp"

namespace jian {

class Lib {
public:
    std::string _lib;

    Lib() : _lib(env("NSP")) {}

    std::string lib() {
        return _lib;
    }

    void lib(const std::string &s) {
        _lib = s;
    }
};

} // namespace jian

#endif

