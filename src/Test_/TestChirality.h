#ifndef TESTCHIRALITY_H
#define TESTCHIRALITY_H

#include "std.h"

class TestChirality {
public:
	TestChirality(char *);
	TestChirality(string);
	void run();

private:
	void init(string);
	string filename;
	double binWidth;
	int maxDist[6];
	int minDist[6];
	int bins;
	Matr_ *chir;
};

#endif
