#include <numeric>

#include "DHMC.hpp"

namespace jian {
	MvEl::MvEl(int a, int b, MvEl::mvel_t t) : type(t) {
		range.push_back({ a, b });
	}

	MvEl::MvEl(int a, int b, int c, int d, MvEl::mvel_t t) : type(t) {
		range.push_back({ a, b });
		range.push_back({ c, d });
	}

	MvEl::MvEl(const helix &h) : type(MvEl::MVEL_HL) {
		int a, b, c, d;

		a = h.head->res1.num - 1;
		d = h.head->res2.num - 1;
		HELIX_EACH(h,
			if (BP->next == NULL) {
				b = BP->res1.num - 1;
				c = BP->res2.num - 1;
			}
		);
		range.push_back({ a, b });
		range.push_back({ c, d });
	}

	MvEl::MvEl(loop *l, MvEl::mvel_t t) : type(t) {
		int a, b, c, d;

		if (t == MVEL_HP) {
			a = l->s.head->res1.num - 1;
			b = l->s.head->res2.num - 1;
			range.push_back({ a, b });
		}
		else if (t == MVEL_IL) {
			a = l->s.head->res1.num - 1;
			b = l->son->s.head->res1.num - 2;
			c = l->son->s.head->res2.num;
			d = l->s.head->res2.num - 1;
			range.push_back({ a, b });
			range.push_back({ c, d });
		}
		else {
			throw "jian::MvEl error!";
		}
	}

	bool MvEl::operator ==(const MvEl &el) const {
		return type == el.type && range == el.range;
	}

	bool MvEl::operator !=(const MvEl &el) const {
		return !(*this == el);
	}

	MvEl *MvEl::operator +(const MvEl &el) const {
		if (el.range.size() == 1) {
			return new MvEl(range[0][0], range[1][1], el.type);
		}
		else if (el.range.size() == 2) {
			return new MvEl(range[0][0], el.range[0][1], el.range[1][0], range[1][1], type);
		}
		else {
			throw "jian::MvEl *jian::MvEl::operator +(const jian::MvEl &el) error!";
		}
	}

	std::ostream &operator <<(std::ostream &stream, const MvEl &el) {
		stream << 
			(el.type == MvEl::MVEL_HL ? "Helix" :
			(el.type == MvEl::MVEL_HP ? "Hairpin" :
			(el.type == MvEl::MVEL_IL ? "Internal Loop" :
			(el.type == MvEl::MVEL_FG ? "Fragment" : "Others")))) << ' ';
		for (auto && frag : el.range) {
			stream << frag[0] << '-' << frag[1] << ' ';
		}
		return stream;
	}

	int MvEl::min() const {
		return std::min_element(range.begin(), range.end(), [](const frag_t &f1, const frag_t &f2) {
			return f1[0] <= f2[0];
		})->at(0);
	}

	int MvEl::max() const {
		return std::max_element(range.begin(), range.end(), [](const frag_t &f1, const frag_t &f2) {
			return f1[1] <= f2[1];
		})->at(1);
	}

	bool MvEl::contains(const MvEl &el) const {
		return std::all_of(el.range.begin(), el.range.end(), [this](const frag_t &f1) {
			return std::any_of(this->range.begin(), this->range.end(), [&f1](const frag_t &f) {
				return f[0] <= f1[0] && f[1] >= f1[1];
			});
		});
	}

	bool MvEl::nips(const MvEl &el) const {
		return range.size() == 2 && range[0][1] + 1 == el.min() && range[1][0] - 1 == el.max();
	}

	void MvEl::merge(std::deque<MvEl *> &dq) {
		LOG << "# Merge ranges..." << std::endl;
		std::deque<MvEl *> els;
		int flag = 1;
		std::map<MvEl *, bool> m;

		while (flag != 0) {
			flag = 0;
			m.clear();
			for (auto && el : dq) m[el] = true;

			auto it = dq.begin();
			for (auto it1 = it; it1 != dq.end(); it1++) {
				if (!m[*it1]) continue;
				for (auto it2 = it1 + 1; it2 != dq.end(); it2++) {
					if (!m[*it2]) continue;
					if ((*it1)->type != MVEL_FG && (*it2)->type != MVEL_FG) {
						MvEl &el1 = *(*it1);
						MvEl &el2 = *(*it2);
						if (el1.contains(el2)) {
							m[*it2] = false;
							flag++;
						}
						else if (el2.contains(el1)) {
							m[*it1] = false;
							flag++;
						}
						else if (el1.nips(el2)) {
							m[*it1] = false;
							m[*it2] = false;
							els.push_back(el1 + el2);
							flag++;
						}
						else if (el2.nips(el1)) {
							m[*it1] = false;
							m[*it2] = false;
							els.push_back(el2 + el1);
							flag++;
						}
					}
				}
			}

			for (auto && el : dq) {
				if (m[el]) {
					els.push_back(el);
				}
				else {
					delete el;
				}
			}
			dq = els;
			els.clear();
		}
	}

	void DHMC::init(const Par &par) {
		nuc3d::mc::MCSM::init(par);

		LOG << "# Set bps" << std::endl;
		set_bps();

		//_ss.resize(_seq.size());
		//std::fill(_ss.begin(), _ss.end(), '.');
		_ss.resize(_seq.size());
		for (auto && c : _ss) c = '.';

		LOG << "# Set 2D trees" << std::endl;
		set_trees();

		LOG << "# Set fragments" << std::endl;
		m_sample_frag = par.has("frag");
		if (m_sample_frag) {
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
		auto & keys = NASS::instance().paired_keys;
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

	void DHMC::set_pseudo_knots_helix(const helix & h) {
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
	}

	void DHMC::set_pseudo_knots() {
		auto it = m_trees.begin();

		// 设置第一个二级结构树中长度为1的helix
		LOOP_TRAVERSE((*it)->head(),
			if (L->has_helix() && L->s.len() == 1) {
				set_pseudo_knots_helix(L->s);
			}
		);

		// 设置除了第一个二级结构树以外的其它的树中的helix
		for (it = m_trees.begin() + 1; it != m_trees.end(); it++) {
			LOOP_TRAVERSE((*it)->head(),
				if (L->has_helix()) {
					set_pseudo_knots_helix(L->s);
				}
			);
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
		for (auto it = m_trees.begin(); it != m_trees.end(); it++) {
			LOOP_TRAVERSE((*it)->head(),
				if (L->has_helix()) {
					restore_helix(&(L->s));
				}
			);
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
		int i;
		m_is_free.resize(_seq.size(), true);

		auto add_el = [this](MvEl *el) {
			for (auto && frag : el->range) {
				for (int i = frag[0]; i <= frag[1]; i++) {
					this->m_is_free[i] = false;
				}
			}
			m_mvels.push_back(el);
		};

		auto set_res_module_types_ss = [&](loop *l, bool is_first) {
			LOOP_TRAVERSE(l,
				if (!_sample_hp && is_first && is_hp(L)) {
					add_el(new MvEl(L, MvEl::MVEL_HP));
				}
				else if (!_sample_il && is_first && is_il(L)) {
					add_el(new MvEl(L, MvEl::MVEL_IL));
				}
				else if (L->has_helix()) {
					add_el(new MvEl(L->s));
				}
			);
		};

		for (auto it = m_trees.begin(); it != m_trees.end(); it++) {
			set_res_module_types_ss((*it)->head(), it == m_trees.begin());
			//std::cout << "m_mvels.size(): " << m_mvels.size() << std::endl;
		}

		std::vector<int> w(m_frag_size);
		for (i = 0; i + m_frag_size - 1 < _seq.size(); i++) {
			std::iota(w.begin(), w.end(), i);
			if (std::any_of(w.begin(), w.end(), [this](int i) {
				return this->m_is_free[i];
			})) {
				m_mvels.push_back(new MvEl(i, i + m_frag_size - 1, MvEl::MVEL_FG));
			}
		}

		//LOG << "information of free:" << std::endl;
		//for (i = 0; i < _seq.size(); i++) {
		//	LOG << i << ' ' << m_is_free[i] << std::endl;
		//}
		//LOG << "before merge: " << std::endl;
		//print_mvels();
		//LOG << "after merge: " << std::endl;

		MvEl::merge(m_mvels);
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
		geom::Superposition sp(m1, m2);
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
		for (auto it = m_trees.begin(); it != m_trees.end(); it++) {
			LOOP_TRAVERSE((*it)->head(),
				if (L->has_helix()) {
					save_helix(&(L->s));
				}
			);
		}
	}

	void DHMC::before_run() {
		LOG << "# Set pseudo-knots" << std::endl;
		set_pseudo_knots();
	}

	void DHMC::mc_select() {
		m_selected_mvel = m_mvels[int(rand() * m_mvels.size())];
	}

	bool DHMC::is_selected(const int &i) const {
		return m_selected_mvel->has(i);
	}

	Vec DHMC::rotating_center() const {
		int beg = m_selected_mvel->min();
		int end = m_selected_mvel->max();
		auto &r1 = _pred_chain[beg];
		auto &r2 = _pred_chain[end];
		Vec origin = Vec::Zero(3);
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
		return origin;
	}


}