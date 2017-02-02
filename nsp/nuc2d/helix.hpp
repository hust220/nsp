#pragma once

#include "jian/pp.hpp"
#include <iostream>
#include <string>
#include <list>
#include <sstream>
#include "bp.hpp"
#include "jian/utils/string.hpp"
#include "jian/utils/ListRange.hpp"

BEGIN_JN

class Helix : public SList<bp> {
public:
	using Self = Helix;
	using Base = SList<bp>;

    Str ss() const {
		STD_ ostringstream stream;
		for (auto &&node : *this) stream << '(';
		for (auto &&node : *this) stream << ')';
		return stream.str();
    }

    Str seq() const {
		int l = 2*size();
		Str str(l, 'X');
		int n = 0;
		for (auto && bp : *this) {
			str[n] = bp.res1.name;
			str[l - 1 - n] = bp.res2.name;
			n++;
		}
        return str;
    }

    Li nums() const {
        Li nums;
		int n = 0;
		for (auto &&bp : *this) {
			nums.insert(std::next(nums.begin(), n), { bp.res1.num - 1, bp.res2.num - 1 });
			n++;
		}
        return nums;
    }

    operator Str() const {
        STD_ ostringstream stream;
        stream << seq() << ' ' << ss();
        for (auto && i : nums()) {
            stream << ' ' << i+1;
        }
        return stream.str();
    }

    friend STD_ ostream &operator <<(STD_ ostream &stream, Helix h) {
        stream << "Helix: " << ' ' << (Str)h;
		return stream;
    }

};

END_JN

