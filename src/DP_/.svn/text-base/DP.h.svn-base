#ifndef DP_H
#define DP_H

#include "../Tools"
class DP {
public:
	DP(string);
	double setEn(int, int);
	double setF(int);
	void run();

protected:
	virtual void setTypes() = 0;
	virtual double getPairEn() = 0;
	double getEn();

	int len_;
	string seq_;
	int *types_;
	int max_span_;
	int min_hairpin_len_;
	double *en_;
	double *f_;
};





#endif




