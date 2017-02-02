#include <jian/utils/exception.hpp>
#include <jntest/UnitTest.hpp>
#include <jnbio/nuc3d/Assemble.hpp>

BEGIN_JN

TEST_CASE(test_assemble)
{
    nuc3d::Assemble ass(Par("seq", "AAAAUUUU")("ss", "(((())))")("name", ""));
    ass.select_templates();
    ass.assemble();
    TEST_CHECK(seq(ass._pred_chain) == "AAAAUUUU");
}

END_JN

