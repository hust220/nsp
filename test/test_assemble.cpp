#include <jntest/UnitTest.hpp>
#include <jnbio/nuc3d/Assemble.hpp>

BEGIN_JN

TEST_CASE(test_assemble)
{
    std::cout << Par("seq", "AAAAUUUU")("ss", "(((())))")("name", "") << std::endl;
    nuc3d::Assemble ass(Par("seq", "AAAAUUUU")("ss", "(((())))")("name", ""));
	ass.predict();
    std::cout << ass._pred_chain << std::endl;
	TEST_CHECK(seq(ass._pred_chain) == "AAAAUUUU");
}

END_JN

