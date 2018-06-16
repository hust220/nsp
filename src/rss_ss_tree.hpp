#pragma once

#include <list>
#include "string.hpp"
#include "rss_sse.hpp"

namespace jian {

class SSTreeRg : 
	public BasicTreeRg<SSTreeRg, PlainTreeIt<PlainTreeEl<SSE>>>
{
public:
	using El = PlainTreeEl<SSE>;
	using Rg = SSTreeRg;
	using Data = SSE;
	using It = PlainTreeIt<El>;
	using Base = BasicTreeRg<Rg, It>;

	Str ss() const {
		Str ss;
		auto p = begin()->head_tail();
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
		auto p = begin()->head_tail();
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

	friend std::ostream &operator <<(std::ostream &stream, SSTreeRg sst) {
		stream << "SST (Secondary Structure Tree) : " << &sst << std::endl;
		for (auto && sse : sst) {
			stream << sse << std::endl;
		}
		return stream;
	}

};

class SSTree : 
	public TreeNt<SSTreeRg> 
{
public:
	SSTree() = default;

	SSTree(const Str &seq, const Str &ss, int hinge = 2);

	// make tree with no broken tag
	void make(const Str &seq, const Str &ss, int hinge = 2);

	// make tree with broken tag
	void make_b(const Str &seq, const Str &ss, int hinge = 2);

};

SSTree ss_tree(Str seq, Str ss, Int h = 1);

}

