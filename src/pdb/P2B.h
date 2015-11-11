#ifndef P2B_H
#define P2B_H

#include "RNA.h"

namespace jian {

class P2B {
public:
	P2B(RNA *);
	void analyze();
	string run();

private:
	int len;
	MatrixXf b;
	MatrixXf n;
	MatrixXf d;
	string ss;
};

} /// namespace jian

#endif

