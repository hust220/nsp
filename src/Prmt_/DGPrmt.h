#ifndef DGPRMT_H
#define DGPRMT_H

#include "std.h"

class DGPrmt {
public:
	DGPrmt(char *);
	DGPrmt(string);
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
