#ifndef JIAN_UTIL_LIB
#define JIAN_UTIL_LIB

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








#endif

