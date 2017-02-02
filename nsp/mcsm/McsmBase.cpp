#include <functional>
#include "../nuc3d/Assemble.hpp"
#include "../nuc3d/Convert.hpp"
#include "McsmBase.hpp"

#define JN_MCXP_TEMPPAR_SET(a) temp_par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
#define JN_MCXP_PAR_SET(a) par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
#define JN_MCXP_TEMP(a) log << PP_STRING3(PP_CAT(mc_, a)) << ' ' << PP_CAT(_mc_, a) << std::endl;

BEGIN_JN

MvEl::MvEl(int a, int b, MvEl::Type t) : type(t) {
	range.push_back({ a, b });
}

MvEl::MvEl(int a, int b, int c, int d, MvEl::Type t) : type(t) {
	range.push_back({ a, b });
	range.push_back({ c, d });
}

MvEl::MvEl(const Helix &h) : type(MvEl::MVEL_HL) {
	int a, b, c, d;

	a = h.front().res1.num - 1;
	d = h.front().res2.num - 1;
	for (auto && bp : h) {
		if (bp.next == NULL) {
			b = bp.res1.num - 1;
			c = bp.res2.num - 1;
		}
	}
	range.push_back({ a, b });
	range.push_back({ c, d });
}

MvEl::MvEl(SSTree::El *l, MvEl::Type t) : type(t) {
	int a, b, c, d;

	if (t == MVEL_HP) {
		a = l->data.helix.front().res1.num - 1;
		b = l->data.helix.front().res2.num - 1;
		range.push_back({ a, b });
	}
	else if (t == MVEL_IL) {
		a = l->data.helix.front().res1.num - 1;
		b = l->son->data.helix.front().res1.num - 2;
		c = l->son->data.helix.front().res2.num;
		d = l->data.helix.front().res2.num - 1;
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
		stream << '(' << frag[0] << '-' << frag[1] << ')';
	}
	return stream;
}

int MvEl::min() const {
	return std::min_element(range.begin(), range.end(), [](const Frag &f1, const Frag &f2) {
		return f1[0] <= f2[0];
	})->at(0);
}

int MvEl::max() const {
	return std::max_element(range.begin(), range.end(), [](const Frag &f1, const Frag &f2) {
		return f1[1] <= f2[1];
	})->at(1);
}

bool MvEl::contains(const MvEl &el) const {
	return std::all_of(el.range.begin(), el.range.end(), [this](const Frag &f1) {
		return std::any_of(this->range.begin(), this->range.end(), [&f1](const Frag &f) {
			return f[0] <= f1[0] && f[1] >= f1[1];
		});
	});
}

bool MvEl::nips(const MvEl &el) const {
	return range.size() == 2 && range[0][1] + 1 == el.min() && range[1][0] - 1 == el.max();
}

void MvEl::merge(Deque<MvEl *> &dq) {
	//log << "# Merge ranges..." << std::endl;
	Deque<MvEl *> els;
	int flag = 1;
	Map<MvEl *, Bool> m;

	if (dq.empty()) return;

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

void MCBase::set_traj_name() {
#ifdef JN_PARA
	m_traj = (mpi_size() == 1 ? to_str(_name, ".traj.pdb") : to_str(_name, ".", mpi_rank() + 1, ".traj.pdb"));
#else
	m_traj = to_str(_name, ".traj.pdb");
#endif
	log << "# Trajectory file: " << m_traj << std::endl;

}

void MCBase::init(const Par &par) {
	TSP::init(par);

	m_selected_mvel = NULL;
	m_sample_mode = SAMPLE_TREE;
	m_cal_en_constraints = true;
	m_max_angle = PI * 0.5;
	m_box = 2;
	m_box_size = 12;
	m_will_write_traj = !_name.empty();

	log << "# Extract residue conformations..." << std::endl;
	for_each_model(to_str(Env::lib(), "/RNA/pars/cg/CG2AA/templates.pdb"), [this](const Model &m, int n) {
		ResConf::extract(m_res_confs, m.residues());
	});
	//MolReader reader(to_str(Env::lib(), "/RNA/pars/cg/CG2AA/templates.pdb"));
	//for (auto it = reader.model_begin(); it != reader.model_end(); it++) {
	//	STD_ cerr << *it << STD_ endl;
	//	ResConf::extract(m_res_confs, it->residues());
	//}
	for (auto && pair : m_res_confs) {
		for (auto && conf : pair.second) {
			conf.res = m_cg->to_cg(conf.res);
		}
	}

	if (m_will_write_traj) {
		set_traj_name();
	}

	log << "# Read parameters..." << std::endl;
	m_par_file = m_cg_type;
	_par->set(m_par_file, "par_file");
	read_parameters();

	log << "# Set parameters..." << std::endl;
	set_parameters(*_par);

	log << "# Print parameters..." << std::endl;
	print_parameters();

	log << "# Set continuous points..." << std::endl;
	set_continuous_pts();

	log << "# Print continuous, angel, dihedral points..." << std::endl;
	for (auto && i : m_continuous_pts) log << i << ' '; log << std::endl;
	for (auto && i : m_ang_pts) log << i << ' '; log << std::endl;
	for (auto && i : m_dih_pts) log << i << ' '; log << std::endl;

	log << "# Read initial structure" << std::endl;
	par.set(m_init_sfile, "init");
	if (m_init_sfile.empty()) {
		if (par.has("init:chain")) {
			_pred_chain = BuildChain()(_seq.size()).m_chain;
		}
		else if (par.has("init:raw")) {
			nuc3d::Assemble assemble(Par(par)("seq", _seq)("ss", _ss)("loop_building", "partial_raw")("name", ""));
			assemble.predict_one();
			_pred_chain = assemble._pred_chain;
		}
	}
	else {
		chain_read_model(_pred_chain, m_init_sfile);
	}

	log << "# Read alignment file" << STD_ endl;
	par.set(m_alignfile, "align");
	if (!m_alignfile.empty()) {
		set_align();
		align_to_fixed_areas();
	}

	log << "# Check constraints" << std::endl;
	validate_constraints();

	//_mc_queue = "heat:30000:20+cool:1000000";
	par.set(_mc_queue, "queue");
}

void MCBase::set_align() {
	Str seq1, seq2;
	Int i, n1, n2, l;
	if (!m_alignfile.empty()) {
		for (auto &&it : FileLines(m_alignfile)) {
			if (it.n == 0) seq1 = it.line;
			else if (it.n == 1) seq2 = it.line;
		}
		if (size(seq1) != size(seq2)) die("Bad align file! The lengths of the two aligned sequences should be equal!");
		l = size(seq1);
		for (i = 0, n1 = 0, n2 = 0; i < l; i++) {
			if (seq1[i] != '-' && seq2[i] != '-') m_align.push_back({ n1, n2 });
			if (seq1[i] != '-') n1++;
			if (seq2[i] != '-') n2++;
		}
	}
#ifndef NDEBUG
	log << "# Print align" << STD_ endl;
	for (auto && pair : m_align) {
		log << pair[0] << ' ' << pair[1] << STD_ endl;
	}
#endif
}

void MCBase::align_to_fixed_areas() {
	Vector<Bool> v(_seq.size(), false);
	for (auto && pair : m_align) {
		v[pair[1]] = true;
	}
	m_fixed_areas.push_back(STD_ move(v));
}

void MCBase::mc_next_step() {
	_mc_step++;
	if (_mc_step >= 100000) {
		//m_sample_mode = SAMPLE_TREE;
		m_cal_en_constraints = true;
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
	log << "Reading parameters of " << m_par_file << "..." << std::endl;
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

	const auto &bks = NASS::instance().break_keys;

    // set broken points
    // from chain lens
	for (n = 0, i = 0; i < m_chain_lens.size(); i++) {
		n += m_chain_lens[i];
		m_brk_pts.push_back(n - 1);
	}
    // from secondary structure
    for (n = 0, i = 0; i < _ss.size(); i++) {
        if (std::find(bks.begin(), bks.end(), _ss[i]) == bks.end()) n++;
        else if (std::find(m_brk_pts.begin(), m_brk_pts.end(), n - 1) == m_brk_pts.end()) m_brk_pts.push_back(n - 1);
    }

    // set continuous points
	for (i = 0; i < _seq.size() - 1; i++) {
		if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end()) {
			m_continuous_pts.push_back(i);
		}
	}
	for (i = 0; i < _seq.size() - 2; i++) {
		if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
			std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 1) == m_brk_pts.end()) {
			m_ang_pts.push_back(i);
		}
	}
	for (i = 0; i < _seq.size() - 3; i++) {
		if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
			std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 1) == m_brk_pts.end() &&
			std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 2) == m_brk_pts.end()) {
			m_dih_pts.push_back(i);
		}
	}
}

void MCBase::mc_write() {
	if (m_will_write_traj) {
		write_traj();
	}
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
		vec_set(axis[0], residue[1]);
		if (t == 1 || t == 3)  vec_set(axis[1], std::plus<>{}, residue[3], residue[5]);
		else if (t == 0 || t == 2) vec_set(axis[1], residue[5]);
		else throw "jian::MCBase::mc_sample_res::get_base_axis error!";
		return axis;
	};

	std::vector<std::function<void()>> actions{
		[this]() {
			int min = m_selected_mvel->min();
			int max = m_selected_mvel->max();
			if (min == 0 && max == _seq.size() - 1) {
				return;
			}
			else if (min == max) {
				ResConf::Confs &confs = m_res_confs[_pred_chain[min].name];
				int l = size(confs);
				int n = int(rand()*l);
				geom::Superposition<Num> sp;
				if (min == 0) {
					Mat x(1, 3);
					mat_set_rows(x, 0, _pred_chain[min + 1]["P"]);
					sp.init(confs[n].p2, x);
				}
				else if (max == _seq.size() - 1) {
					Mat x(1, 3);
					mat_set_rows(x, 0, _pred_chain[min]["P"]);
					sp.init(confs[n].p1, x);
				}
				else {
					Mat x(2, 3);
					mat_set_rows(x, 0, _pred_chain[min]["P"], _pred_chain[min + 1]["P"]);
					sp.init(confs[n].pp, x);
				}
				Residue r = confs[n].res;
				for (auto && atom : r) sp.apply(atom);
				set_atoms(_pred_chain[min], r);
			}
			else {
				Num d1 = -1;
				Num d2 = -1;
				if (min > 0) d1 = geom::distance(_pred_chain[min - 1][1], _pred_chain[min][0]);
				if (max < size(_seq) - 1) d2 = geom::distance(_pred_chain[max][1], _pred_chain[max + 1][0]);
				/*if (d1 > 4.0 || d2 > 4.0) {
					int index = int(rand() * 3);
					Num dist = (rand() - 0.5) * 0.3 * _mc_max_shift;
					for (int i = min; i <= max; i++) {
						for (auto && atom : _pred_chain[i]) {
							atom[index] += dist;
						}
						space_update_item(i);
					}
				}
				else */if (min == 0 || max == size(_seq) - 1 || d1 > 4.0 || d2 > 4.0) {
					int t = ((min == 0 || d1 > 4.0) ? max : min);
					int index = int(rand() * 3);
					double dih = (rand() - 0.5) * m_max_angle;
					auto &&rot = geom::rot_mat(index, dih);
					Vec origin(3);
					vec_set(origin, _pred_chain[t][0]);
					for (int i = min; i <= max; i++) {
						for (auto && atom : _pred_chain[i]) {
							geom::rotate(atom, origin, rot);
						}
						space_update_item(i);
					}
				}
				else {
					geom::RotateAlong<double> rotate_along(_pred_chain[min][0], _pred_chain[max + 1][0], m_max_angle * (rand() - 0.5));
					for (int i = min; i <= max; i++) {
						for (auto && atom : _pred_chain[i]) {
							rotate_along(atom);
						}
						space_update_item(i);
					}
				}
			}
		},

		// rotate base
		[this, &get_base_axis]() {
			for (int i = 0; i < _seq.size(); i++) {
				if (is_selected(i)) {
					auto axis = get_base_axis(_pred_chain[i]);
					geom::RotateAlong<double> rotate_along(axis[0], axis[1], m_max_angle * (rand() - 0.5));
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
			double dih = (rand() - 0.5) * m_max_angle;
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
	actions[0]();

	//int min = m_selected_mvel->min();
	//int max = m_selected_mvel->max();
	//if (m_sample_mode == SAMPLE_SSE) {
	//	if (min == max) {
	//		actions[std::vector<int>{0, 1, 2}[int(rand() * 3)]]();
	//	}
	//	else {
	//		actions[std::vector<int>{0, 2}[int(rand()*2)]]();
	//	}
	//}
	//else if (m_sample_mode == SAMPLE_TREE) {
	//	if (min == max) {
	//		actions[std::vector<int>{0, 1}[int(rand() * 2)]]();
	//	}
	//	else {
	//		actions[0]();
	//	}
	//}
	//else {
	//	throw "jian::MCBase::sample_res error!";
	//}

}

void MCBase::mc_back() {
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

void MCBase::backup() {
	_moved_atoms.clear();

		for (int i = 0; i < _seq.size(); i++) {
			if (is_selected(i)) {
				for (auto && atom : _pred_chain[i]) {
					_moved_atoms.push_back(atom);
				}
			}
		}
}

void MCBase::init_space() {
	m_item_space.resize(_seq.size());
	for (int i = 0; i < _seq.size(); i++) {
		SpaceVal &s = space_val(i);
		s.push_back(i);
		m_item_space[i] = &s;
	}
}

int MCBase::space_index(double n) const {
	return int((n + 1000) / m_box_size);
}

MCBase::SpaceItem &MCBase::item(int i) {
	return _pred_chain[i][2];
}

MCBase::SpaceVal &MCBase::space_val(int i) {
	SpaceItem &a = item(i);
	//std::cout << a << std::endl;
	//std::cout << space_index(a[0]) << ' ' << space_index(a[1]) << ' ' << space_index(a[2]) << std::endl;
	return m_space[space_index(a[0])][space_index(a[1])][space_index(a[2])];
}

void MCBase::space_update_item(int i) {
	SpaceVal &n = space_val(i);
	SpaceVal &o = *(m_item_space[i]);
	if (&o != &n) {
		o.erase(std::find(o.begin(), o.end(), i));
		m_item_space[i] = &n;
		n.push_back(i);
	}
}

void MCBase::restore_raw() {
	Mat dist;
	Int l = size(_seq);
	dg_dist_init(dist, l);
	//dg_dist_read_ss(dist, _seq, _ss);
	Vector<Int> v(l);
	for (Int i = 0; i < l; i++) {
		auto it = STD_ find_if(m_align.begin(), m_align.end(), [i](auto &&p) {return p[0] == i; });
		if (it == m_align.end()) v[i] = -1;
		else v[i] = it->at(1);
	}
	dg_dist_read_chain(dist, _pred_chain, v);
	Mat && c = DG(dist)();
	auto cg = CG::fac_t::make("1p");
	_pred_chain = cg->to_aa(c, 0, c.rows() - 1);
}

void MCBase::run() {
	log << "# Check initial structure..." << std::endl;
	if (_pred_chain.empty()) die("Please give me an initial structure for optimizing!");
	if (!m_alignfile.empty()) restore_raw();

	log << "# Initializing running..." << std::endl;
	before_run();

	save_fixed_ranges();

	log << "# Carrying on CG processing with the Chain..." << std::endl;
	_pred_chain = m_cg->to_cg(_pred_chain);

	log << "# Init space..." << std::endl;
	init_space();

	log << "# MC..." << std::endl;
	mc_run();

	log << "# Print Constraints and Distances..." << std::endl;
	print_final_constraints();

	log << "# Finishing running..." << std::endl;
	finish_run();

	log << "# Coarsed Grained To All Atom..." << std::endl;
	_pred_chain = m_cg->to_aa(_pred_chain);

	log << "# Restore helix..." << std::endl;
	restore_fixed_ranges();

	log << "# Transform..." << std::endl;
	transform();

}

void MCBase::transform() {
	int i = 0;
	for (Residue &r : _pred_chain) {
		r = convert_res(r, to_str(_seq[i]));
		i++;
	}
}

void MCBase::print_final_constraints() {
	double d;
	int i, j;
	log << "Print Contacts:" << std::endl;
	for (auto && c : _constraints.contacts) {
		i = c.key[0];
		j = c.key[1];
		d = dist_two_res(_pred_chain[i], _pred_chain[j]);
		log << i << ' ' << j << " weight:" << c.weight << " dist:" << d << std::endl;
	}
	log << "Print Distances:" << std::endl;
	for (auto && c : _constraints.distances) {
		i = c.key[0];
		j = c.key[1];
		d = dist_two_res(_pred_chain[i], _pred_chain[j]);
		log << i << ' ' << j << " value:" << c.value << " weight:" << c.weight << " dist:" << d << std::endl;
	}
}

void MCBase::before_run() {}
void MCBase::finish_run() {}

Str MCBase::file_parameters() const {
	return "3drna";
}

void MCBase::save_fixed_ranges() {}
void MCBase::restore_fixed_ranges() {}

END_JN
