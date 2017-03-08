#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../mc.hpp"
#include "../pdb.hpp"
#include "jian/geom.hpp"
#include "../nuc2d.hpp"
#include "../cg.hpp"
#include "jian/pp.hpp"
#include "../mcsm.hpp"
#include "../cg/ResFrags.hpp"
#include "jian/utils/file.hpp"
#include "../nuc3d/BuildHelix.hpp"
#include "../nuc3d/transform.hpp"
#include "../nuc3d/TemplRec.hpp"
#include "../nuc3d/TSP.hpp"

BEGIN_JN

class DHMC : public MCSM {
public:
	using Range = Ai<4>;

	Deque<SP<SSTree>> m_trees;
	Map<SSE *, Chain *> m_saved_helices;
	Int m_frag_size;
	ResFrags *m_frags;
    Bool m_save_bg;
	Bool m_sample_frag;
	Bool m_sample_all_res;
	Bool m_not_sample_hp;
	Bool m_not_sample_il;
	Bool m_set_mvel_pk;
	Bool m_all_free;
	Map<Str, Array<Num, 6>> m_bp_distances;
//	Vector<Bool> m_is_free;
//    Vector<Bool> m_is_single_pair;

	DHMC() = default;

	~DHMC();

	static SP<DHMC> make(const Par &par) {
		SP<DHMC> dhmc = STD_ make_shared<DHMC>();
		dhmc->init(par);
		return dhmc;
	}

	void init(const Par &par);

	void bps_to_constraints();

	void set_base_mvels();

	void remove_useless_constraints();

	void print_constraints();

	void print_mvels();

	void set_trees();

	void set_bps();

	void set_pseudo_knots();

	void transform_saved_helix(Chain *h, const Li & nums);

	void restore_helix(SSE *sse);

	virtual void restore_fixed_ranges();

	// MC related methods

	virtual void mc_sample();

	void save_helix(SSE *sse);

	virtual void save_fixed_ranges();

	virtual void before_run();

	virtual void mc_select();

	virtual bool is_selected(const I &i) const;

	virtual Vec rotating_center() const;

};

//template<typename T>
//void chain_refine(Chain &chain, SSTree::El *l, const Frags &fixed_ranges = {}, S traj = "") {
//	Par par;
//	par._orig_pars = { "nsp", "" };
//	Str seq = SSTree::make_range(l).seq();
//	Str ss = SSTree::make_range(l).ss();
//	par._pars["traj"].push_front(traj);
//	par._pars["seq"].push_front(seq);
//	par._pars["ss"].push_front(ss);
//	DHMC mc;
//	mc.init(par);
//	mc._pred_chain = chain;
//	mc.run();
//	chain = mc._pred_chain;
//}
//
END_JN

