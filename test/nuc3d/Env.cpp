#include "UnitTest.hpp"
#include <jian/nuc3d/Env.hpp>
#include <iostream>
#include <cstring>

TEST_CASE(Env_lib) {
    using namespace jian;
    TEST_CHECK(!Env::lib().empty());
}


