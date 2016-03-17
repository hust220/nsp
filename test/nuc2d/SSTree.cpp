#include "UnitTest.hpp"
#include <jian/nuc2d/SSTree.hpp>
#include <iostream>
#include <cstring>

TEST_CASE(SSTree_make) {
    using namespace jian;
    SSTree ss_tree;
    ss_tree.make("AAAACCCCUUUU", "((((....))))");
    LOOP_TRAVERSE(ss_tree.head,
        L->print();
    );
//    TEST_CHECK(i == 71);
}


