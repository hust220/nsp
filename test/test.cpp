#include <iostream>
#include <jian/pp/map.hpp>

int main(int argc, char **argv) {
#define TEST(a) int a
#define LS 1, 2, 3
	std::cout << PP_STRING3(JN_MAP(TEST, LS)) << std::endl;
#define AA(...) __VA_ARGS__
	//std::cout << PP_STRING3() << std::endl;
	return 0;
}