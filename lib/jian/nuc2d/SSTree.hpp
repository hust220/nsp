#pragma once

#include <list>
#include "../utils/string.hpp"
#include "SSE.hpp"

BEGIN_JN

class SSTreeRange : 
	public TreeRange<SSE>
{
public:
	using Base = TreeRange<SSE>;

	SSTreeRange(SSE *node = NULL) : Base(node) {}

	Str ss() const {
		Str ss;
		auto p = head->head_tail();
		ss.resize(p.second - p.first + 1);
		for (auto &&sse : *this) {
			for (auto && bp : sse.helix) {
				ss[bp.res1.num - 1] = bp.res1.type;
				ss[bp.res2.num - 1] = bp.res2.type;
			}
			for (auto && res : sse.loop) {
				ss[res.num - 1] = res.type;
			}
		}
		return ss;
	}

	Str seq() const {
		Str seq;
		auto p = head->head_tail();
		seq.resize(p.second - p.first + 1);
		for (auto &&sse : *this) {
			for (auto &&bp : sse.helix) {
				seq[bp.res1.num - 1] = bp.res1.name;
				seq[bp.res2.num - 1] = bp.res2.name;
			}
			for (auto &&res : sse.loop) {
				seq[res.num - 1] = res.name;
			}
		}
		return seq;
	}

	friend STD_ ostream &operator <<(STD_ ostream &stream, SSTreeRange sst) {
		for (auto && sse : sst) {
			stream << sse << STD_ endl;
		}
	}

};

class SSTree : public Tree<SSTreeRange> {
public:
	SSTree() = default;

	SSTree(const Str &seq, const Str &ss, int hinge = 2);

	// make tree with no broken tag
	void make(const Str &seq, const Str &ss, int hinge = 2);

	// make tree with broken tag
	void make_b(const Str &seq, const Str &ss, int hinge = 2);

};

END_JN

