#ifndef STPRMT2_H
#define STPRMT2_H

#include "std.h"

class StPrmt2 {
public:
	StPrmt2(char *);
	StPrmt2(string);
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
