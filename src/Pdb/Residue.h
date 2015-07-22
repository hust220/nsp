#ifndef RESIDUE_H_INCLUDED
#define RESIDUE_H_INCLUDED

#include "Atom.h"

namespace jian {

class Residue
{
public:
	Residue() {}
	Residue(Residue *r) : num(r->num), number(r->number), atomNum(r->atomNum), 
	                      rnaName(r->rnaName), name(r->name), atoms(r->atoms) {}
	Residue(const Residue &r) : num(r.num), number(r.number), atomNum(r.atomNum), 
	                            rnaName(r.rnaName), name(r.name), atoms(r.atoms) {}
	Residue &operator =(const Residue &r) {
		num = r.num;
		number = r.number;
		atomNum = r.atomNum;
		rnaName = r.rnaName;
		name = r.name;
		atoms = r.atoms;
		return *this;
	}
	friend ostream &operator <<(ostream &, const Residue &);
	Residue(vector<string> &, string);
	Residue(Point *, int, int);
	Atom &operator [](int);
	const Atom &operator [](int) const;
    Atom &operator [](std::string);
    const Atom &operator [](std::string) const;
    std::vector<Atom> operator [](const std::vector<std::string> &) const;
	Point *getBaseMassCenter();
	int getAtomPos(string, Point &);
	int getAmount();
	Point *getBaseVec();
	int nextTo(Residue &);

	vector<Atom>::iterator begin();
	vector<Atom>::iterator end();

	int num;
	string number;
	int atomNum;
	string rnaName;
	string name = "X";
	vector<Atom> atoms;
};

} /// namespace jian

#endif // RESIDUE_H_INCLUDED

