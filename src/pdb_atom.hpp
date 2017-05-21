#pragma once

#include <string>
#include <array>
#include <map>
#include "matrix.hpp"

BEGIN_JN

class Atom : public std::array<double, 3> {
public:
	S name;
	int num;
	double mass;

	JN_DEFAULT_CONSTRUCTORS(Atom);

	Atom(S name, double x, double y, double z, int num = -1) {
		init(name, x, y, z, num);
	}

	void init(S name, double x, double y, double z, int num = -1);

	void set_name(const S &s);
	void set_mass();

};

#define JN_DEF_ATOMS\
	refs<Atom> atoms() { return refs<Atom>().append(*this); }\
	refs<const Atom> atoms() const { return refs<const Atom>().append(*this); }

END_JN

