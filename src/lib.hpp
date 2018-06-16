#ifndef JIAN_UTIL_LIB
#define JIAN_UTIL_LIB

#include <string>

namespace jian {

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

}

#endif

