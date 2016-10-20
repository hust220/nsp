#pragma once

#include <string>
#include <exception>

namespace jian {

class Error : public std::exception {
private:
    std::string _inf;

public:
    Error(const std::string &inf = "") : _inf(inf) {}

    virtual const char *what() const noexcept {
        return _inf.c_str();
    }
};

}

