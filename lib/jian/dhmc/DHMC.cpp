#include <numeric>
#include "../nuc2d/NASS.hpp"
#include "DHMC.hpp"

namespace jian {
	void DHMC::init(const Par &par) {
		MCSM::init(par);

		m_set_mvel_pk = false;
		m_pk_ahead = par.has("pk_ahead");
		m_sample_frag = par.has("frag");
		m_sample_all_res = par.has("sample_all_res");
		m_not_sample_hp = par.has("not_sample_hp");
		m_not_sample_il = par.has("not_sample_il");
		m_all_free = par.has("all_free");

		LOG << "# Set bps" << std::endl;
		set_bps();

		LOG << "# Transform bps to constraints..." << std::endl;
		bps_to_constraints();

		//if (m_pk_ahead) {
		//	// no operation
		//}
		//else if (m_sample_all_res) {
		//	for (auto && c : _ss) {
		//		c = '.';
		//	}
		//}
		//else {
		//	for (auto && c : _ss) {
		//		if (c != '.' && c != '(' && c != ')') {
		//			c = '.';
		//		}
		//	}
		//}

		LOG << "# Set 2D trees" << std::endl;
		set_trees();

		if (m_sample_frag) {
			LOG << "# Set fragments" << std::endl;
			m_frag_size = 3;
			par.set(m_frag_size, "frag_size");
			m_frags = &(Frags::instance(m_cg->m_cg, m_frag_size));
		}
		else {
			m_frag_size = 1;
		}

		LOG << "# Set moving elements" << std::endl;
		set_mvels();

		LOG << "# Set moving elements of each base" << std::endl;
		set_base_mvels();

		LOG << "# Print moving elements" << std::endl;
		print_mvels();

		//throw "debugging...";

		LOG << "# Remove useless constraints" << std::endl;
		remove_useless_constraints();

		LOG << "# Print constraints" << std::endl;
		print_constraints();

	}

	DHMC::~DHMC() {
		for (auto && el : m_mvels) {
			delete el;
		}
		for (auto && h : m_saved_helices) {
			delete h.second;
		}
	}

	void DHMC::set_base_mvels() {
		m_base_mvels.resize(_seq.size());
		for (auto && el : m_mvels) {
			for (auto && frag : el->range) {
				for (int i = frag[0]; i <= frag[1]; i++) {
					m_base_mvels[i] = el;
				}
			}
		}
	}

	void DHMC::remove_useless_constraints() {
		Constraints cs;
		for (auto && c : _constraints.contacts) {
			if (m_base_mvels[c.key[0]] != m_base_mvels[c.key[1]]) {
				cs.contacts.push_back(c);
			}
		}
		for (auto && c : _constraints.distances) {
			if (m_base_mvels[c.key[0]] != m_base_mvels[c.key[1]]) {
				cs.distances.push_back(c);
			}
		}
		_constraints = cs;
	}

	void DHMC::print_constraints() {
		for (auto && c : _constraints.distances) {
			LOG << c.key[0] << ' ' << c.key[1] << std::endl;
		}
	}

	void DHMC::print_mvels() {
		for (auto && el : m_mvels) {
			LOG << *el << std::endl;
		}
	}

	void DHMC::set_trees() {
		m_trees.push_back(std::make_shared<SSTree>());
		m_trees.back()->make_b(_seq, _ss, 2);
		const NASS::PairedKeys & keys = NASS::instance().paired_keys;
		for (auto it = keys.begin() + 1; it != keys.end(); it++) {
			auto ss = partial_ss(_ss, *it);
			if (std::any_of(ss.begin(), ss.end(), [](auto &&c) {return c != '.' && c != '&'; })) {
				m_trees.push_back(std::make_shared<SSTree>());
				m_trees.back()->make_b(_seq, ss, 1);
			}
			else {
				break;
			}
		}
	}

	void DHMC::set_bps() {
		m_bps = NASS::get_bps(_ss);
	}

	void DHMC::bps_to_constraints() {
		int i, j;
		std::vector<val_t> dists{ 18.2, 15.2, 10.7, 4.2, 6.3, 6.4 };

		auto foo = [this, &dists](int n1, int n2) {
			//int t1 = pdb::res_type(_pred_chain[n1].name);
			//int t2 = pdb::res_type(_pred_chain[n2].name);

			for (int i = 0; i < 6; i++) {
				m_distance_constraints.push_back({ { n1, i },{ n2, i }, dists[i], dists[i] });
			}

			//int a = ((t1 == 0 || t1 == 2) ? 5 : 4);
			//int b = ((t2 == 0 || t2 == 2) ? 5 : 4);
			//m_distance_constraints.push_back({ { n1, 3 },{ n2, 3 }, 3.9, 4.5 });
			//m_distance_constraints.push_back({ { n1, a },{ n2, b }, 3.9, 4.5 });
		};

		for (i = 0; i < _seq.size(); i++) {
			j = m_bps[i];
			if (i < j) {
				foo(i, j);
			}
		}
	}

	void DHMC::translate_pseudo_knots_helix(Model & m, const std::list<int> & nums) {
		int n = m_cg->res_size() * nums.size() / 2;
		Mat x(n, 3), y(n, 3);
		int i = 0;
		int l = 0;
		for (auto && j : nums) {
			if (i < nums.size() / 2) {
				Residue && res1 = m_cg->to_cg(m[0][i]);
				Residue && res2 = m_cg->to_cg(_pred_chain[j]);
				for (int k = 0; k < m_cg->res_size(); k++) {
					for (int t = 0; t < 3; t++) {
						x(l, t) = res1[k][t];
						y(l, t) = res2[k][t];
					}
					l++;
				}
			}
			i++;
		}
		auto sp = geom::suppos(x, y);
		INIT_SUPPOS(sp);
		for (auto && res : m[0]) {
			for (auto && atom : res) {
				APPLY_SUPPOS(atom, sp);
			}
		}
	}

	void DHMC::set_pseudo_knots() {
		if (m_pk_ahead) {
			LOG << "# Set pseudo-knots" << std::endl;

			auto it = m_trees.begin();

			auto foo = [this](const helix &h) {
				auto && seq = h.seq();
				auto && m = build_helix(seq);
				auto && nums = h.nums();

				assert(nums.size() >= 2 && nums.size() % 2 == 0);
				translate_pseudo_knots_helix(m, nums);

				int i = 0;
				for (auto && n : nums) {
					_pred_chain[n] = std::move(m[0][i]);
					i++;
				}
			};

			// 设置第一个二级结构树中长度为1的helix
			LOOP_TRAVERSE((*it)->head(),
				if (L->has_helix() && L->s.len() == 1) {
					foo(L->s);
				}
			);

			// 设置除了第一个二级结构树以外的其它的树中的helix
			for (it = m_trees.begin() + 1; it != m_trees.end(); it++) {
				LOOP_TRAVERSE((*it)->head(),
					if (L->has_helix()) {
						foo(L->s);
					}
				);
			}
		}
	}

	void DHMC::transform_saved_helix(Chain *h, const std::list<int> & nums) {
		int n = m_cg->res_size() * nums.size();
		Mat x(n, 3), y(n, 3);
		int i = 0;
		int l = 0;
		for (auto && j : nums) {
			Residue && res1 = m_cg->to_cg(h->at(i));
			Residue && res2 = m_cg->to_cg(_pred_chain[j]);
			for (int k = 0; k < m_cg->res_size(); k++) {
				for (int t = 0; t < 3; t++) {
					x(l, t) = res1[k][t];
					y(l, t) = res2[k][t];
				}
				l++;
			}
			i++;
		}
		auto sp = geom::suppos(x, y);
		INIT_SUPPOS(sp);
		for (auto && res : *h) {
			for (auto && atom : res) {
				APPLY_SUPPOS(atom, sp);
			}
		}
	}

	void DHMC::restore_helix(helix *h) {
		auto && seq = h->seq();
		auto && nums = h->nums();
		Chain *saved = m_saved_helices[h];

		transform_saved_helix(saved, nums);

		int i = 0;
		for (auto && n : nums) {
			_pred_chain[n] = std::move(saved->at(i));
			i++;
		}
	}

	void DHMC::restore_fixed_ranges() {
		if (!m_all_free) {
			LOG << "# " << __FUNCTION__ << std::endl;

			auto end = (m_pk_ahead ? m_trees.end() : std::next(m_trees.begin()));
			for (auto it = m_trees.begin(); it != end; it++) {
				LOOP_TRAVERSE((*it)->head(),
					if (L->has_helix()) {
						restore_helix(&(L->s));
					}
				);
			}
		}
	}

	bool DHMC::is_hp(loop *l) {
		if (l->has_loop() && l->has_helix() && l->num_sons() == 0) {
			int flag = 0, n = 0;
			LOOP_EACH(l,
				char &c = _ss[RES->num - 1];
				if (c != '.'  && c != '(' && c != ')') {
					flag = 1;
					break;
				}
				else {
					n++;
				}
			);
			return flag == 0 && n <= 14;
		}
		else {
			return false;
		}
	};

	bool DHMC::is_il(loop *l) {
		if (l->has_loop() && l->has_helix() && l->num_sons() == 1) {
			int flag = 0, n = 0;
			LOOP_EACH(l,
				char &c = _ss[RES->num - 1];
				if (c != '.'  && c != '(' && c != ')') {
					flag = 1;
					break;
				}
				else {
					n++;
				}
			);
			return flag == 0 && n <= 20;
		}
		else {
			return false;
		}
	};

	void DHMC::set_mvels() {
		//int i;
		//m_is_free.resize(_seq.size(), true);

		//auto add_el = [this](MvEl *el) {
		//	for (auto && frag : el->range) {
		//		for (int i = frag[0]; i <= frag[1]; i++) {
		//			this->m_is_free[i] = false;
		//		}
		//	}
		//	m_mvels.push_back(el);
		//};

		//auto set_res_module_types_ss = [&](loop *l, bool is_first) {
		//	LOOP_TRAVERSE(l,
		//		if (m_not_sample_hp && is_first && is_hp(L)) {
		//			add_el(new MvEl(L, MvEl::MVEL_HP));
		//		}
		//		else if (m_not_sample_il && is_first && is_il(L)) {
		//			add_el(new MvEl(L, MvEl::MVEL_IL));
		//		}
		//		else if (L->has_helix()) {
		//			add_el(new MvEl(L->s));
		//		}
		//	);
		//};

		//if (m_sample_all_res) {
		//	// pass
		//}
		//else if (m_set_mvel_pk) {
		//	for (auto it = m_trees.begin(); it != m_trees.end(); it++) {
		//		set_res_module_types_ss((*it)->head(), it == m_trees.begin());
		//	}
		//}
		//else {
		//	set_res_module_types_ss(m_trees.front()->head(), true);
		//}

		//MvEl::merge(m_mvels);

		//int max_extend_len = 5;
		//for (auto && el : m_mvels) {
		//	int min = el->min();
		//	int max = el->max();
		//	for (int i = 1; i <= max_extend_len; i++) {
		//		if (min - i >= 0 && m_is_free[min - i]) {
		//			m_mvels.push_back(new MvEl(min - i, max, MvEl::MVEL_FG));
		//		}
		//		else {
		//			break;
		//		}
		//	}
		//	for (int i = 1; i <= max_extend_len; i++) {
		//		if (max + i < _seq.size() && m_is_free[max + i]) {
		//			m_mvels.push_back(new MvEl(min, max + i, MvEl::MVEL_FG));
		//		}
		//		else {
		//			break;
		//		}
		//	}
		//}

		//for (int frag_size = 1; frag_size <= 10; frag_size++) {
		//	std::vector<int> w(frag_size);
		//	for (i = 0; i + frag_size - 1 < _seq.size(); i++) {
		//		std::iota(w.begin(), w.end(), i);
		//		if (std::all_of(w.begin(), w.end(), [this](int i) {
		//			return this->m_is_free[i];
		//		})) {
		//			m_mvels.push_back(new MvEl(i, i + frag_size - 1, MvEl::MVEL_FG));
		//		}
		//	}
		//}

		const NASS::PairedKeys &keys = NASS::instance().paired_keys;
		auto is_paired = [this](int i, int j) {return m_bps[i] == j && (_ss[i] == '(' || m_pk_ahead); };
		auto foo = [this, &is_paired, &keys](int i, int j) {
			if (is_paired(i, j) && i > 0 && j < size(_seq) - 1 && is_paired(i-1, j+1)) return false;
			auto end = (m_pk_ahead ? keys.end() : std::next(keys.begin(), 1));
			std::vector<int> v(std::distance(keys.begin(), end));
			//int s = 0;
			int m;
			for (int n = i; n <= j; n++) {
				char c = _ss[n];
				auto it1 = std::find_if(keys.begin(), end, [&c](const NASS::PairedKey &pair) {
					return pair.first == c;
				});
				auto it2 = std::find_if(keys.begin(), end, [&c](const NASS::PairedKey &pair) {
					return pair.second == c;
				});
				if (it1 != end) {
					m = std::distance(keys.begin(), it1);
					v[m]++;
					//s++;
				}
				else if (it2 != end) {
					m = std::distance(keys.begin(), it2);
					v[m]--;
					//s--;
					if (v[m] < 0) return false;
				}
			}
			return std::all_of(v.begin(), v.end(), [](int s) {return s == 0; });
			//return s == 0;
		};

		for (int i = 0; i < _seq.size(); i++) {
			for (int j = i; j < _seq.size(); j++) {
				if (m_all_free || foo(i, j)) {
					m_mvels.push_back(new MvEl(i, j, MvEl::MVEL_FG));
				}
			}
		}

	}

	// MC related methods

	void DHMC::mc_sample() {
		if (m_selected_mvel->type == MvEl::MVEL_FG && m_sample_frag) {
			//LOG << "sample frag" << std::endl;
			mc_sample_frag();
		}
		else {
			//LOG << "sample res" << std::endl;
			mc_sample_res();
		}
	}

	void DHMC::mc_sample_frag() {
		int i, j, k, id;
		Mat m, m1, m2;
		Mat *mat;

		if (m_selected_mvel->type != MvEl::MVEL_FG)
			throw "jian::DHMC::mc_sample_frag error!";

		backup();

		std::ostringstream stream;
		MvEl::frag_t &f = m_selected_mvel->range[0];
		for (i = f[0]; i <= f[1]; i++) {
			stream << _seq[i];
		}

		auto &ids = m_frags->m_ids[stream.str()];
		id = ids[int(ids.size() * rand())];
		mat = m_frags->m_mats[id];
		//Chain *c1 = m_frags->m_chains_aa[id];
		//Chain *c2 = m_frags->m_chains_cg[id];
		//static int n = 0;
		//if (n < 10) {
		//	mol_write(*c1, "aa.frag.aa." + JN_STR(n + 1) + ".pdb");
		//	mol_write(*c2, "aa.frag.cg." + JN_STR(n + 1) + ".pdb");
		//	std::ostringstream str;
		//	str << "aa.mat." << n + 1 << ".txt";
		//	std::ofstream ofile(str.str().c_str());
		//	ofile << *m << std::endl;
		//	ofile.close();
		//}
		//n++;

		m = *mat;
		m1.resize(2 * m_cg->res_size(), 3);
		m2.resize(2 * m_cg->res_size(), 3);
		//Chain c;
			//c.push_back(_pred_chain[f[0] + i]);
		for (j = 0; j < m_cg->res_size(); j++) {
			for (k = 0; k < 3; k++) {
				m1(j, k) = m(j, k);
				m1(m_cg->res_size() + j, k) = m((m_frag_size - 1) * m_cg->res_size() + j, k);
				m2(j, k) = _pred_chain[f[0]][j][k];
				m2(m_cg->res_size() + j, k) = _pred_chain[f[0] + (m_frag_size - 1)][j][k];
			}
		}
		//static int n = 0;
		//if (n < 10) mol_write(c, "aa.frag." + JN_STR(n + 1) + ".pdb");
		//n++;
		geom::Superposition<double> sp(m1, m2);
		//std::cout << "before sp: " << sp.rmsd << std::endl;
		sp.apply_m(m);
		//std::cout << "after sp: " << geom::rmsd(m1, m2) << std::endl;
		for (i = 0; i < m_frag_size; i++) {
			if (m_is_free[f[0] + i]) {
				for (j = 0; j < m_cg->res_size(); j++) {
					for (k = 0; k < 3; k++) {
						_pred_chain[f[0] + i][j][k] = m(i * m_cg->res_size() + j, k);
					}
				}
				space_update_item(f[0] + i);
			}
		}
	}

	void DHMC::save_helix(helix *h) {
		auto && nums = h->nums();

		Chain *c = new Chain;

		for (auto && n : nums) {
			c->push_back(_pred_chain[n]);
		}

		m_saved_helices[h] = c;
	}

	void DHMC::save_fixed_ranges() {
		if (!m_all_free) {
			LOG << "# " << __FUNCTION__ << std::endl;

			auto end = (m_pk_ahead ? m_trees.end() : std::next(m_trees.begin()));
			for (auto it = m_trees.begin(); it != end; it++) {
				LOOP_TRAVERSE((*it)->head(),
					if (L->has_helix()) {
						save_helix(&(L->s));
					}
				);
			}
		}
	}

	void DHMC::before_run() {
		set_pseudo_knots();
	}

	void DHMC::mc_select() {
		m_selected_mvel = m_mvels[int(rand() * m_mvels.size())];
	}

	bool DHMC::is_selected(const int &i) const {
		if (m_selected_mvel == NULL) {
			return false;
		}
		else if (m_sample_mode == SAMPLE_SSE) {
			return m_selected_mvel->has(i);
		}
		else if (m_sample_mode == SAMPLE_TREE) {
			return m_selected_mvel->minmax_has(i);
		}
		else {
			throw "jian::DHMC::is_selected error! Illegal sample mode!";
		}
	}

	Vec DHMC::rotating_center() const {
		Vec origin = Vec::Zero(3);

		if (m_selected_mvel->type == MvEl::MVEL_FG) {
			for (int i = 0; i < 3; i++) origin[i] = _pred_chain[m_selected_mvel->range[0][0]][0][i];
		}
		else {
			int beg = m_selected_mvel->min();
			int end = m_selected_mvel->max();
			auto &r1 = _pred_chain[beg];
			auto &r2 = _pred_chain[end];
			int n_atom = 0;
			for (auto && atom : r1) {
				for (int i = 0; i < 3; i++) origin[i] += atom[i];
				n_atom++;
			}
			for (auto && atom : r2) {
				for (int i = 0; i < 3; i++) origin[i] += atom[i];
				n_atom++;
			}
			for (int i = 0; i < 3; i++) origin[i] /= n_atom;
		}
		return origin;
	}


}