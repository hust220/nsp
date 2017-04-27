#ifndef JIAN_UTIL_LIB
#define JIAN_UTIL_LIB

#include <string>
#include "../lib.hpp"

BEGIN_JN

class Lib {
public:
    S _lib;

    Lib() : _lib(env("NSP")) {}

    S lib() {
        return _lib;
    }

    void lib(const S &s) {
        _lib = s;
    }
};

END_JN

#endif

