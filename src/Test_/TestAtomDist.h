#ifndef TESTATOMDIST_H
#define TESTATOMDIST_H

#include "std.h"

class TestAtomDist {
public:
	TestAtomDist(char *, char *, int, char *);
	TestAtomDist(string, string, int, string);
	void run();

private:
	void init(string, string, int, string);
	string name1;
	string name2;
	int num;
	string filename;
	double binWidth;
	int maxDist;
	int minDist;
	int maxDist2;
	int minDist2;
	int bins;
	int *distribution;
	int *distribution2;
};

#endif
