#ifndef DISTANAL_H
#define DISTANAL_H

#include "../Pdb.h"

namespace jian {

class DistAnal {
public:
	DistAnal(int = 20, double = 0.5, string = "average");
	~DistAnal();

	void readRNA(const RNA &);
	void train();
	double scoring();

	void readObsParm(string = "");
	void readRefParm(string = "");
	void initObsProb();
	void initRefProb();

	void printObsParm();
	void printRefParm();
	void printObsProb();
	void printRefProb();

	double *getScore();
	double score[5];
	double * nuc_score_;
	int *nuc_len_;
	double * stacking_score_;
	int *stacking_len_;
	Obj<Matr_> pairing_score_;
	Obj<Matr_> pairing_len_;

private:
	string reference_state;
	int *num;
	int *type;
	int *ntLen;
	Point **list;
	int len;

	/* 1: i-i+1; 2: i-i+2; 3: i-i+3; 4: i-i+n */
	int *obs_parm_[5];
	int *ref_parm_[5];
	double *obs_prob_[5];
	double *ref_prob_[5];

	double interval;
	int cutoff;
	int bins;
	double penalty;
};

} /// namespace jian

#endif

