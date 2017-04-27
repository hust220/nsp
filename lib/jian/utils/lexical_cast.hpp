#pragma once

#include <string>
#include <sstream>
#include "traits.hpp"

BEGIN_JN

template<typename T, typename U>
T lexical_cast(U && u) {
    std::stringstream stream;
    T t;

    stream << u;
    stream >> t;
    return t;
}

END_JN

