#include "McsmBase.hpp"

#define JN_MCXP_TEMPPAR_SET(a) temp_par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
#define JN_MCXP_PAR_SET(a) par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
#define JN_MCXP_TEMP(a) LOG << PP_STRING3(PP_CAT(mc_, a)) << ' ' << PP_CAT(_mc_, a) << std::endl;

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



	void MCBase::init(const Par &par) {
		TSP::init(par);

		m_selected_mvel = NULL;
		m_sample_mode = SAMPLE_SSE;

		LOG << "# Set the file of trajectory..." << std::endl;
		std::ostringstream stream;
#ifdef JN_PARA
		stream << _name << "." << g_mpi->m_rank + 1 << ".traj.pdb";
#else
		stream << _name << ".traj.pdb";
#endif
		m_traj = stream.str();

		LOG << "# Set continuous points..." << std::endl;
		set_continuous_pts();

		LOG << "# Print continuous, angel, dihedral points..." << std::endl;
		for (auto && i : m_continuous_pts) LOG << i << ' '; LOG << std::endl;
		for (auto && i : m_ang_pts) LOG << i << ' '; LOG << std::endl;
		for (auto && i : m_dih_pts) LOG << i << ' '; LOG << std::endl;

		LOG << "# Read initial structure" << std::endl;
		if (par.has("pdb")) chain_read_model(_pred_chain, par.get("pdb"));

		LOG << "# Check constraints" << std::endl;
		validate_constraints();

		//_mc_queue = "heat:30000:20+cool:1000000";
		//par.set(_mc_queue, "mc_queue");
	}

	void MCBase::mc_next_step() {
		_mc_step++;
		if (_mc_step >= 100000) {
			m_sample_mode = SAMPLE_TREE;
		}
	}

	void MCBase::validate_constraints() {
		int l = _seq.size();
		int i, j;
		for (auto && ct : _constraints.distances) {
			i = ct.key[0];
			j = ct.key[1];
			if (i < 0 || i >= l || j < 0 || j >= l) {
				std::cerr << "Illegal constraints pair (" << i << ", " << j << ") in your constraints file!" << std::endl;
				std::exit(1);
			}
		}
	}

	void MCBase::read_parameters() {
		LOG << "Reading parameters of " << m_par_file << "..." << std::endl;
		Par temp_par(Env::lib() + "/RNA/pars/nuc3d/mc/" + m_par_file + ".par");
		JN_MAP(JN_MCXP_TEMPPAR_SET, JN_MCXP_PARS1, JN_MCXP_PARS2)
	}

	void MCBase::set_parameters(const Par &par) {
		JN_MAP(JN_MCXP_PAR_SET, JN_MCXP_PARS1, JN_MCXP_PARS2)
	}

	void MCBase::print_parameters() {
		JN_MAP(JN_MCXP_TEMP, JN_MCXP_PARS1, JN_MCXP_PARS2)
	}

	void MCBase::set_continuous_pts() {
		int i, n;

		n = 0;
		for (i = 0; i < m_chain_lens.size(); i++) {
			n += m_chain_lens[i];
			m_brk_pts.push_back(n - 1);
		}
		for (int i = 0; i < _seq.size() - 1; i++) {
			if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end()) {
				m_continuous_pts.push_back(i);
			}
		}
		for (int i = 0; i < _seq.size() - 2; i++) {
			if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
				std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 1) == m_brk_pts.end()) {
				m_ang_pts.push_back(i);
			}
		}
		for (int i = 0; i < _seq.size() - 3; i++) {
			if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
				std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 1) == m_brk_pts.end() &&
				std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 2) == m_brk_pts.end()) {
				m_dih_pts.push_back(i);
			}
		}
	}

	void MCBase::mc_write() {
		write_traj();
		write_en();
		mc_num_writing++;
	}

	void MCBase::write_traj() {
		if (mc_num_writing == 1) {
			file::clean(m_traj);
		}
		std::ofstream output(m_traj.c_str(), std::ios::app);
		m_writer.bind_stream(output);
		m_writer.write_model([&]() {
			this->m_writer.write(this->_pred_chain);
		});
		output.close();
	}

	void MCBase::mc_sample() {
		mc_sample_res();
	}

	void MCBase::mc_sample_res() {
		auto get_base_axis = [](const Residue &residue) {
			int t = pdb::res_type(residue.name);
			std::array<Vec, 2> axis;
			axis[0].resize(3);
			axis[1].resize(3);
			for (int i = 0; i < 3; i++) {
				axis[0][i] = residue[1][i];
			}
			if (t == 1 || t == 3) {
				for (int i = 0; i < 3; i++) {
					axis[1][i] = residue[3][i] + residue[5][i];
				}
			}
			else if (t == 0 || t == 2) {
				for (int i = 0; i < 3; i++) {
					axis[1][i] = residue[5][i];
				}
			}
			else {
				throw "jian::MCBase::mc_sample_res::get_base_axis error!";
			}
			return axis;
		};

		std::vector<std::function<void()>> actions{
			// rotate along P-P
			[this]() {
				int min = m_selected_mvel->min();
				int max = m_selected_mvel->max();
				if (max + 1 < _seq.size()) {
					geom::RotateAlong<double> rotate_along(_pred_chain[min][0], _pred_chain[max+1][0], PI * 0.1 * (rand() - 0.5));
					for (int i = min; i <= max; i++) {
						for (auto && atom : _pred_chain[i]) {
							rotate_along(atom);
						}
						space_update_item(i);
					}
				}
				else {
					int index = int(rand() * 3);
					double dih = (rand() - 0.5) * PI * 0.1;
					auto &&rot = geom::rot_mat(index, dih);
					Vec origin(3);
					for (int i = 0; i < 3; i++) origin[i] = _pred_chain[min][0][i];
					for (int i = min; i <= max; i++) {
						for (auto && atom : _pred_chain[i]) {
							geom::rotate(atom, origin, rot);
						}
						space_update_item(i);
					}
				}
			},

			// rotate base
			[this, &get_base_axis]() {
				for (int i = 0; i < _seq.size(); i++) {
					if (is_selected(i)) {
						auto axis = get_base_axis(_pred_chain[i]);
						geom::RotateAlong<double> rotate_along(axis[0], axis[1], PI * 0.1 * (rand() - 0.5));
						for (int j = 3; j < 6; j++) {
							rotate_along(_pred_chain[i][j]);
						}
						space_update_item(i);
					}
				}
			},
				// translate
			[this]() {
				int index = int(rand() * 3);
				double dist = (rand() - 0.5) * 1 * _mc_max_shift;
				for (int i = 0; i < _seq.size(); i++) {
					if (is_selected(i)) {
						for (auto && atom : _pred_chain[i]) {
							atom[index] += dist;
						}
						space_update_item(i);
					}
				}
			},

			// rotate fixed in P
			[this]() {
				int index = int(rand() * 3);
				double dih = (rand() - 0.5) * PI * 0.1;
				auto &&rot = geom::rot_mat(index, dih);
				auto &&origin = rotating_center();
				for (int i = 0; i < _seq.size(); i++) {
					if (is_selected(i)) {
						for (auto && atom : _pred_chain[i]) {
							geom::rotate(atom, origin, rot);
						}
						space_update_item(i);
					}
				}
			}

		};

		backup();

		int min = m_selected_mvel->min();
		int max = m_selected_mvel->max();
		if (m_sample_mode == SAMPLE_SSE) {
			if (min == max) {
				actions[2+int(2 * rand())]();
			}
			else {
				actions[1+int(3 * rand())]();
			}
		}
		else if (m_sample_mode == SAMPLE_TREE) {
			if (min == max) {
				actions[(rand() < 0.5 ? 0 : 1)]();
			}
			else {
				actions[0]();
			}
		}
		else {
			throw "jian::MCBase::sample_res error!";
		}

	}

	void MCBase::mc_back() {
		if (m_sample_mode == SAMPLE_TREE) {
			for (int i = m_selected_mvel->min(); i <= m_selected_mvel->max(); i++) {
				for (auto && atom : _pred_chain[i]) {
					atom = _moved_atoms.front();
					_moved_atoms.pop_front();
				}
				space_update_item(i);
			}
		}
		else if (m_sample_mode == SAMPLE_SSE) {
			for (int i = 0; i < _seq.size(); i++) {
				if (is_selected(i)) {
					for (auto && atom : _pred_chain[i]) {
						atom = _moved_atoms.front();
						_moved_atoms.pop_front();
					}
					space_update_item(i);
				}
			}
		}
		else {
			throw "jian::MCBase::mc_back error!";
		}
	}

	void MCBase::backup() {
		_moved_atoms.clear();

		if (m_sample_mode == SAMPLE_TREE) {
			for (int i = m_selected_mvel->min(); i <= m_selected_mvel->max(); i++) {
				for (auto && atom : _pred_chain[i]) {
					_moved_atoms.push_back(atom);
				}
			}
		}
		else if (m_sample_mode == SAMPLE_SSE) {
			for (int i = 0; i < _seq.size(); i++) {
				if (is_selected(i)) {
					for (auto && atom : _pred_chain[i]) {
						_moved_atoms.push_back(atom);
					}
				}
			}
		}
		else {
			throw "jian::MCBase::backup error!";
		}
	}

	void MCBase::init_space() {
		m_item_space.resize(_seq.size());
		for (int i = 0; i < _seq.size(); i++) {
			space_val_t &s = space_val(i);
			s.push_back(i);
			m_item_space[i] = &s;
		}
	}

	int MCBase::space_index(double n) const {
		return int((n + 1000) / m_box_size);
	}

	MCBase::item_t &MCBase::item(int i) {
		return _pred_chain[i][0];
	}

	MCBase::space_val_t &MCBase::space_val(int i) {
		item_t &a = item(i);
		return m_space[space_index(a[0])][space_index(a[1])][space_index(a[2])];
	}

	void MCBase::space_update_item(int i) {
		space_val_t &n = space_val(i);
		space_val_t &o = *(m_item_space[i]);
		if (&o != &n) {
			o.erase(std::find(o.begin(), o.end(), i));
			m_item_space[i] = &n;
			n.push_back(i);
		}
	}

	void MCBase::run() {
		LOG << "# Check initial structure..." << std::endl;
		if (num_residues(_pred_chain) == 0) {
			throw "Please give an initial structure before the optimization procedure!";
		}

		LOG << "# Read parameters..." << std::endl;
		m_par_file = m_cg_type;
		_par->set(m_par_file, "par_file");
		read_parameters();

		LOG << "# Set parameters..." << std::endl;
		set_parameters(*_par);

		LOG << "# Print parameters..." << std::endl;
		print_parameters();

		LOG << "# Initializing running..." << std::endl;
		before_run();

		LOG << "# Save helices" << std::endl;
		save_fixed_ranges();

		LOG << "# Carrying on CG processing with the Chain..." << std::endl;
		_pred_chain = m_cg->to_cg(_pred_chain);

		LOG << "# Init space..." << std::endl;
		init_space();

		LOG << "# MC..." << std::endl;
		mc_run();

		LOG << "# Print Constraints and Distances..." << std::endl;
		print_final_constraints();

		LOG << "# Finishing running..." << std::endl;
		finish_run();

		LOG << "# Coarsed Grained To All Atom..." << std::endl;
		cg_to_aa(chain_to_coords());

		LOG << "# Restore helix..." << std::endl;
		restore_fixed_ranges();

		LOG << "# Transform..." << std::endl;
		this->transform();

	}

	void MCBase::transform() {
		Model m;
		m.push_back(_pred_chain);
		_pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
	}

	Mat MCBase::chain_to_coords() {
		int n_atoms = num_atoms(_pred_chain);
		Mat c(n_atoms, 3);
		int n_atom = 0;
		for (int i = 0; i < _pred_chain.size(); i++) {
			for (auto && atom : _pred_chain[i]) {
				for (int j = 0; j < 3; j++) {
					c(n_atom, j) = atom[j];
				}
				n_atom++;
			}
		}
		return c;
	}

	void MCBase::print_final_constraints() {
		double d;
		int i, j;
		LOG << "Print Contacts:" << std::endl;
		for (auto && c : _constraints.contacts) {
			i = c.key[0];
			j = c.key[1];
			d = dist_two_res(_pred_chain[i], _pred_chain[j]);
			LOG << i << ' ' << j << " weight:" << c.weight << " dist:" << d << std::endl;
		}
		LOG << "Print Distances:" << std::endl;
		for (auto && c : _constraints.distances) {
			i = c.key[0];
			j = c.key[1];
			d = dist_two_res(_pred_chain[i], _pred_chain[j]);
			LOG << i << ' ' << j << " value:" << c.value << " weight:" << c.weight << " dist:" << d << std::endl;
		}
	}

	void MCBase::cg_to_aa(const Mat &c) {
		_pred_chain = m_cg->to_aa(c, 0, c.rows() - 1);
	}

	void MCBase::before_run() {}
	void MCBase::finish_run() {}

	std::string MCBase::file_parameters() const {
		return "3drna";
	}

	void MCBase::save_fixed_ranges() {}
	void MCBase::restore_fixed_ranges() {}

}
