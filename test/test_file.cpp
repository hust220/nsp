#include "UnitTest.hpp"
#include <jian/utils/file.hpp>

BEGIN_JN

TEST_CASE(file_test)
{
	STD_ ofstream ofile("file-test.txt");
	ofile << "hi";
	ofile << "hello";
	ofile.close();
	for (auto &&it : JN_ FileLines("file-test.txt")) {
		TEST_CHECK(it.line == "hi");
	}
}

END_JN


