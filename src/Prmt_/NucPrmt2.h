#ifndef NUCPRMT2_H
#define NUCPRMT2_H

#include "std.h"

class NucPrmt2 {
public:
	NucPrmt2(char *);
	NucPrmt2(string);
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
