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

	// MvEl: Moving element
	class MvEl {
	public:
		enum mvel_t {
			MVEL_HL, // helix
			MVEL_HP, // hairpin
			MVEL_IL, // internal loop
			MVEL_FG // fragment
		};
		using frag_t = std::array<int, 2>;
		using range_t = std::vector<frag_t>;

		mvel_t type;
		range_t range;

		MvEl(int a, int b, mvel_t t);

		MvEl(int a, int b, int c, int d, mvel_t t);

		MvEl(const helix &h);

		MvEl(loop *l, mvel_t t);

		int min() const;

		int max() const;

		bool operator ==(const MvEl &el) const;

		bool operator !=(const MvEl &el) const;

		MvEl *operator +(const MvEl &el) const;

		friend std::ostream &operator <<(std::ostream &, const MvEl &el);

		bool contains(const MvEl &el) const;

		bool nips(const MvEl &el) const;

		bool has(int n) const {
			return std::find_if(range.begin(), range.end(), [&n](const frag_t &frag) {
				return frag[0] <= n && n <= frag[1];
			}) != range.end();
		}

		static void merge(std::deque<MvEl *> &dq);
	};

	using fixed_ranges_t = std::list<std::array<int, 2>>;
	using fixed_mvels_t = std::list<MvEl>;

	class DHMC : public nuc3d::mc::MCSM {
	public:
		using range_t = std::array<int, 4>;

		std::deque<std::shared_ptr<SSTree>> m_trees;
		std::map<helix *, Chain *> m_saved_helices;
		std::deque<MvEl *> m_mvels;
		std::vector<MvEl *> m_base_mvels;
		MvEl *m_selected_mvel;
		fixed_mvels_t m_fixed_mvels;
		int m_frag_size;
		Frags *m_frags;
		bool m_sample_frag;
		std::vector<bool> m_is_free;

		template<typename T>
		std::string partial_ss(std::string ss, T &&pair) {
			for (auto && c : ss) {
				c = (c == pair.first ? '(' : (c == pair.second ? ')' : (c == '&' ? '.' : c)));
			}
			return ss;
		}

		DHMC() = default;

		~DHMC();

		void init(const Par &par);

		bool is_hp(loop *l);

		bool is_il(loop *l);

		void set_base_mvels();

		void remove_useless_constraints();

		void print_constraints();

		void print_mvels();

		void set_trees();

		void set_bps();

		void translate_pseudo_knots_helix(Model & m, const std::list<int> & nums);

		void set_pseudo_knots_helix(const helix & h);

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

