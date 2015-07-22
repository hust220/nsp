#ifndef DNA_H
#define DNA_H

#include "Model.h"

namespace jian {


class DNA :public Model {
public:
	DNA() : Model() {}
	DNA(const Model &model): Model(model) {}
	DNA(DNA *dna) {
		name = dna->name;
		len = dna->len;
		chains = dna->chains;
	}
	DNA(const DNA &dna) {
		name = dna.name;
		len = dna.len;
		chains = dna.chains;
	}
	DNA &operator =(const DNA &dna) {
		name = dna.name;
		len = dna.len;
		chains = dna.chains;
		return *this;
	}
	friend ostream &operator <<(ostream &, const DNA &);
	void print();
	void write(string);
	DNA(char *);
	DNA(string);
	DNA *copy();
	void readPDB(string);
	void push(Chain *);
	void push(Chain &);
	Chain &operator [](int);
	void updateChains(string);

	void setLen();
	void setResNum();

	int getLen();
	int totalAtoms();
	double getDist(int, int);
	string getSeq();
	string getChain();
	
	/* IO function */
//	void print();
//	void printAsDNA();
//	void write(string);
	
//	/* assemble function */
//	void move(double, double, double);
//	void rotate(Matr_ *m);
//	void format();
//	void addP();
//	void mutate(string);
//	void rotateByX(double);
//	void rotateByZ(double);

	/* attributes */
	string name;
	int len;
};

} /// namespace jian

//class DNAs {
//public:
//	DNAs();
//	~DNAs();
//
//	int getLen();
//	void resize(int);
//	DNA *at(int);
//	DNA *operator [](int);
//	void push(DNA *);
//	
//	vector<DNA *> DNAList;
//};

#endif // DNA_H

