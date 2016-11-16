#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../mc.hpp"
#include "../utils/file.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "../cg.hpp"
#include "../pp.hpp"
#include "../nuc3d/BuildHelix.hpp"
#include "../nuc3d/transform.hpp"
#include "../nuc3d/TemplRec.hpp"
#include "../nuc3d/TSP.hpp"
#include "../mcsm.hpp"
#include "../cg/Frags.hpp"

namespace jian {

	class DHMC : public MCSM {
	public:
		using range_t = std::array<int, 4>;

		std::deque<std::shared_ptr<SSTree>> m_trees;
		std::map<helix *, Chain *> m_saved_helices;
		int m_frag_size;
		Frags *m_frags;
		bool m_sample_frag;
		bool m_pk_ahead;
		bool m_sample_all_res;
		bool m_not_sample_hp = false;
		bool m_not_sample_il = false;
		bool m_set_mvel_pk = false;
		std::vector<bool> m_is_free;

		template<typename T>
		std::string partial_ss(std::string ss, T &&pair) {
			for (auto && c : ss) {
				c = (c == pair.first ? '(' : (c == pair.second ? ')' : '.'));
				//c = (c == pair.first ? '(' : (c == pair.second ? ')' : (c == '&' ? '.' : c)));
			}
			return ss;
		}

		DHMC() = default;

		~DHMC();

		void init(const Par &par);

		void bps_to_constraints();

		bool is_hp(loop *l);

		bool is_il(loop *l);

		void set_base_mvels();

		void remove_useless_constraints();

		void print_constraints();

		void print_mvels();

		void set_trees();

		void set_bps();

		void translate_pseudo_knots_helix(Model & m, const std::list<int> & nums);

		//void set_pseudo_knots_helix(const helix & h);

		void set_pseudo_knots();

		void transform_saved_helix(Chain *h, const std::list<int> & nums);

		void restore_helix(helix *h);

		virtual void restore_fixed_ranges();

		void set_mvels();

		//bool in_fixed_mvels(const MvEl &el);

		// MC related methods

		virtual void mc_sample();

		void mc_sample_frag();

		void save_helix(helix *h);

		virtual void save_fixed_ranges();

		virtual void before_run();

		virtual void mc_select();

		virtual bool is_selected(const int &i) const;

		virtual Vec rotating_center() const;

	};

	template<typename T>
	void chain_refine(Chain &chain, loop *l, const fixed_ranges_t &fixed_ranges = {}, std::string traj = "") {
		Par par;
		par._orig_pars = { "nsp", "" };
		std::string seq, ss;
		seq_read_tree(seq, l);
		ss_read_tree(ss, l);
		par._pars["traj"].push_front(traj);
		par._pars["seq"].push_front(seq);
		par._pars["ss"].push_front(ss);
		DHMC mc;
		mc.init(par);
		mc._pred_chain = chain;
		//mc.m_fixed_ranges = fixed_ranges;
		mc.run();
		chain = mc._pred_chain;
	}

} // namespace jian

