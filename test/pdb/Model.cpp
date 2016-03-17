#include "UnitTest.hpp"
#include <jian/pdb/Model.hpp>
#include <iostream>
#include <cstring>

TEST_CASE(Model) {
    jian::Model model("/home/jian/1Y26.pdb");
    TEST_CHECK(jian::num_residues(model) == 71);
}


