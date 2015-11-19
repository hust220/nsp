#ifndef NUCPRMT3_H
#define NUCPRMT3_H

#include "std.h"

class NucPrmt3 {
public:
	NucPrmt3(char *);
	NucPrmt3(string);
	void run();

private:
	void init(string);
	string filename;
	int *lens;
	double binWidth;
	int bins;
	Matr_ **prmt;
	Matr_ **prmtMinMax;
};

#endif
