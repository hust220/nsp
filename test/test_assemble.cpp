#include "UnitTest.hpp"
#include "jian/nuc3d/Assemble.hpp"

BEGIN_JN

TEST_CASE(test_assemble)
{
	nuc3d::Assemble ass(Par("seq", "AAAAUUUU")("ss", "(((())))")("name", "aa"));
	ass.predict();
	Model model;
	mol_read(model, "aa.pred.pdb");
	TEST_CHECK(seq(model) == "AAAAUUUU");
}

END_JN


