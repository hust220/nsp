#include "UnitTest.hpp"
#include "jian/nuc3d/Assemble.hpp"

BEGIN_JN

TEST_CASE(test_assemble)
{
    std::cout << Par("seq", "AAAAUUUU")("ss", "(((())))")("name", "") << std::endl;
    nuc3d::Assemble ass(Par("seq", "AAAAUUUU")("ss", "(((())))")("name", ""));
    std::cout << "hi" << std::endl;
	ass.predict();
    std::cout << "hi" << std::endl;
	//Model model;
	//mol_read(model, "aa.pred.pdb");
    std::cout << seq(ass._pred_chain) << std::endl;
	TEST_CHECK(seq(ass._pred_chain) == "AAAAUUUU");
}

END_JN


