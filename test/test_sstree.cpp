#include <jntest/UnitTest.hpp>
#include <jnbio/nuc2d/SSTree.hpp>

BEGIN_JN

TEST_CASE(test_sstree)
{
	SSTree sst("AAAACCCCUUUUAAUU", "((((....))))(())", 2);
	SSTree sst2(sst);
	sst.make("AACCAACCUUCCAAAACCUUCCUUAAUUUUCC", "((..((..))..((((..))..))(())))..", 1);
	TEST_CHECK(size(sst) == 6);
	TEST_CHECK(size(sst2) == 3);
	TEST_CHECK(size(sst.range()) == 6);
}

END_JN


