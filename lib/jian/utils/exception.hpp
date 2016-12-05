#pragma once

#include <string>
#include <exception>

BEGIN_JN

class Error : public std::exception {
private:
    S _inf;

public:
    Error(const S &inf = "") : _inf(inf) {}

    virtual const char *what() const noexcept {
        return _inf.c_str();
    }
};

}

