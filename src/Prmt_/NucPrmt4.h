#ifndef NUCPRMT4_H
#define NUCPRMT4_H

#include "std.h"

class NucPrmt4 {
public:
	NucPrmt4(char *);
	NucPrmt4(string);
	void run();

private:
	void init(string);
	string filename;
	double binWidth;
	int bins;
	Matr_ *prmt;
	Matr_ *prmtMinMax;
};

#endif
