#include <jntest/UnitTest.hpp>
#include <jian/utils/ListRange.hpp>

BEGIN_JN

TEST_CASE(list_test)
{
	SList<int> ls;
	ls.push_back(1, 2, 3, 4, 5);
	Vi v{ 1, 2, 3, 4, 5 };
	Int n = 0;
	for (auto && i : ls) {
		TEST_CHECK(v[n] == i);
		n++;
	}
	n = 0;
	for (auto it = ls.begin(); it != ls.end(); it++) {
		TEST_CHECK(v[n] == *it);
		n++;
	}
	TEST_CHECK(*(ls.begin()) == 1);
	TEST_CHECK(size(ls) == 5);

	SList<Pair<int, Str>> ls2;
	ls2.push_back(STD_ make_pair(1, "a"));
	TEST_CHECK((*ls2.begin()).first == 1);
	TEST_CHECK(ls2.begin()->first == 1);
	TEST_CHECK(ls2.begin()->second == "a");
	TEST_CHECK(ls2.begin()->second == "a");
}

END_JN


