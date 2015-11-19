#ifndef STRAND_H
#define STRAND_H

#include "res.h"
#include "../Pdb.h"

class strand {
public:
	strand();
	double getDist(RNA *);
	int getLen();
	res *head, *tail;
};

class strandInfo {
public:
		string name;
		int len;
		double dist;
		string seq;
		string src;
};

class strandInfoEx {
public:
	strandInfoEx(); 
	strandInfoEx(strandInfo *si, double score); 

	strandInfo *si;
	strandInfoEx *next;
	double score;
};

class strandInfoList {
public:
	strandInfoList(); 
	int getLen();
	void add(strandInfo *hi, double score); 

	strandInfoEx *head;
};

#endif

