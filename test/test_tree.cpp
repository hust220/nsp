#include <jntest/UnitTest.hpp>
#include <jian/utils/Tree.hpp>

BEGIN_JN

TEST_CASE(test_tree)
{
	Tree<int> tree;
	auto root = tree.set_root(1);
	auto son = root->set_son(2);
	auto brother = root->set_brother(3);
	son->set_brother(4);
	brother->set_son(5);
	TEST_CHECK(tree.front() == 1);
	TEST_CHECK(size(tree) == 5);
	Vi v{ 1, 2, 4, 3, 5 };
	int n = 0;
	for (auto &&i : tree) {
		TEST_CHECK(v[n] == i);
		n++;
	}
	for (auto && it : tree.path()) {
		//STD_ cout << it->data << STD_ endl;
	}
}

END_JN


