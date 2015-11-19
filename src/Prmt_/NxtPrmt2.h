#ifndef NXTPRMT2_H
#define NXTPRMT2_H

#include "std.h"

class NxtPrmt2 {
public:
	NxtPrmt2(char *);
	NxtPrmt2(string);
	void run();

private:
	void init(string);
	string filename;
	double binWidth;
	int bins;
	Matr_ **prmt;
	Matr_ **prmtMinMax;
};

#endif
