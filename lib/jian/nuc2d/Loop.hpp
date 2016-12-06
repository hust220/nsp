#pragma once

#include "../utils/traits.hpp"
#include "../utils/ListRange.hpp"
#include "res.hpp"

BEGIN_JN

class LoopRange;
using Loop = ListEntity<LoopRange>;

class LoopRange :
	public ListRange<res>
{
public:
	using Self = LoopRange;
	using Base = ListRange<res>;
	using Node = Base::Node;
	using It = Base::It;

	LoopRange(Node *node = NULL) : Base(node) {}

	Str ss() const {
		STD_ ostringstream stream;
		for (auto &&res : *this) stream << res.type;
		return stream.str();
	}

	Str seq() const {
		STD_ ostringstream stream;
		for (auto &&res : *this) {
			if (res.type != '&') {
				stream << res.type;
			}
		}
		return stream.str();
	}

	operator Str() const {
		std::ostringstream stream;
		stream << seq() << ' ' << ss() << ' ';
		for (auto && res : *this) stream << res.num << ' ';
		return stream.str();
	}

	friend STD_ ostream &operator <<(STD_ ostream &stream, Loop l) {
		stream << "Loop: " << ' ' << (Str)l;
	}

};

END_JN