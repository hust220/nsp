#ifndef POTPRMT_H
#define POTPRMT_H

#include "std.h"

class PotPrmt {
public:
	PotPrmt(char *);
	PotPrmt(string);
	void run();

private:
	void init(string);
	void analyze(Point *, Point *, int, int, int, int);
	string filename;
	double binWidth;
	int bins;
	double **theta;
};

#endif
