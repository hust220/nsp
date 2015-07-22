#ifndef BPPrmt_H
#define BPPrmt_H

#include "std.h"

class BPPrmt {
public:
	BPPrmt(char *);
	BPPrmt(string);
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
