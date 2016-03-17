#include "UnitTest.hpp"
#include <jian/etl/std.hpp>
#include <iostream>
#include <cstring>

TEST_CASE(std_FOR) {
    int n = 0; FOR((i, 3), n += i);
    TEST_CHECK(n == 3);
}

TEST_CASE(std_FOLD) {
    TEST_CHECK(FOLD(_1 + _2, 0, (1, 2, 3)) == 6);
}

