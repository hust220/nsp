#include "ParBp.hpp"
#include "../utils/Par.hpp"
#include "../utils/log.hpp"
#include "Score.hpp"

namespace jian {

	void Score::init() {
		int i, j;

		parbp.reset(new ParBp);
		m_bin = 0.3;
		m_cutoff = 20;
		m_bins = int(std::ceil(m_cutoff / m_bin));
		m_res_size = m_cg->res_size();
		//m_res_size = 6;
		m_num_types = m_res_size * 4;
		m_rows = m_num_types * m_num_types;
		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++) {
				m_counts_bp[i][j] = Mat3i::Zero(40, 40, 8);
				m_freqs_bp[i][j] = Mat3::Zero(40, 40, 8);
				m_weights_bp[i][j] = 0;
				m_counts_st53[i][j] = Mat3i::Zero(20, 20, 20);
				m_freqs_st53[i][j] = Mat3::Zero(20, 20, 20);
				m_weights_st53[i][j] = 0;
				m_counts_st35[i][j] = Mat3i::Zero(20, 20, 20);
				m_freqs_st35[i][j] = Mat3::Zero(20, 20, 20);
				m_weights_st35[i][j] = 0;
			}
		}
		read_counts();
		read_pars();
		set_freqs();
	}

	void Score::run(const Chain & c) {
		int i, j, l;
		Chain chain;

		chain = m_cg->to_cg(c);
		l = chain.size();
		m_score = 0;
		for (i = 0; i < l; i++) {
			for (j = i + 1; j < l; j++) {
				m_score += en_pairing(chain[i], chain[j]);
			}
		}
	}

	void Score::train(const Chain &c) {
		int i, j, len;
		Chain chain;
		std::vector<int> t;

		auto counts_bp_add = [this](int i, int j, auto && v) {
			if (std::fabs(v[0]) < 10 && std::fabs(v[1]) < 10 && std::fabs(v[2]) < 2) {
				int a = int((v[0] + 10) / 0.5);
				int b = int((v[1] + 10) / 0.5);
				int c = int((v[2] + 2) / 0.5);
				this->m_counts_bp[i][j](a, b, c)++;
			}
		};

		auto counts_st_add = [this](auto && m, int i, int j, auto && v) {
			if (std::fabs(v[0]) < 5 && std::fabs(v[1]) < 5 && std::fabs(v[2]) < 5) {
				int a = int((v[0] + 5) / 0.5);
				int b = int((v[1] + 5) / 0.5);
				int c = int((v[2] + 5) / 0.5);
				m[i][j](a, b, c)++;
			}
		};

		chain = m_cg->to_cg(c);
		m_chain = &chain;
		set_indices();
		len = chain.size();

		t.resize(len);
		for (i = 0; i < len; i++) {
			t[i] = pdb::res_type(chain[i].name);
			if (t[i] > 3 || t[i] < 0) throw "jian::Score::train error!";
		}
		for (i = 0; i < len; i++) {
			for (j = i + 1; j < len; j++) {
				parbp->anal(chain[i], chain[j]);
				if (parbp->is_paired()) {
					counts_bp_add(t[i], t[j], parbp->o21_);
					counts_bp_add(t[j], t[i], parbp->o12_);
				}
				else if (parbp->is_stacked()) {
					counts_st_add(m_counts_st53, t[i], t[j], parbp->o21_);
					counts_st_add(m_counts_st35, t[j], t[i], parbp->o12_);
				}
			}
		}

		for (i = 0; i < len - 1; i++) {
			update_counts_stacking(i, i + 1);
			for (j = i + 2; j < len; j++) {
				update_counts_pairing(i, j);
			}
		}
	}

	void Score::set_indices() {
		int i, len;

		len = m_chain->size();
		m_indices.resize(len);
		for (i = 0; i < len; i++) {
			m_indices[i] = m_map[m_chain->at(i).name];
		}
	}

	void Score::update_counts_stacking(int n1, int n2) {
		int i, j;
		double d;
		int n, a, b;

		a = m_indices[n1];
		b = m_indices[n2];
		const Residue &r1 = m_chain->at(n1);
		const Residue &r2 = m_chain->at(n2);
		ParBp parbp(r1, r2);
		if (parbp.is_stacked()) {
			for (i = 0; i < m_res_size; i++) {
				for (j = 0; j < m_res_size; j++) {
					d = geom::distance(r1[i], r2[j]);
					if (d < 20) {
						n = (a * m_res_size + i) * m_num_types + (b * m_res_size + j);
						m_counts_stacking(n, int(d / m_bin))++;
					}
				}
			}
		}
	}

	void Score::update_counts_pairing(int n1, int n2) {
		int i, j;
		double d;
		int n, a = m_indices[n1], b = m_indices[n2];
		const Residue &r1 = m_chain->at(n1);
		const Residue &r2 = m_chain->at(n2);
		ParBp parbp(r1, r2);
		bool is_wc = parbp.is_wc();
		bool is_nwc = parbp.is_nwc();
		for (i = 0; i < m_res_size; i++) {
			for (j = 0; j < m_res_size; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (d < 20) {
					n = (a * m_res_size + i) * m_num_types + (b * m_res_size + j);
					m_counts_pairing(n, int(d / m_bin))++;
					if (is_wc) {
						m_counts_wc(n, int(d / m_bin))++;
					}
					if (is_nwc) {
						m_counts_nwc(n, int(d / m_bin))++;
					}
				}
			}
		}
	}

	void Score::read_counts() {
		std::ifstream ifile;

		auto foo = [this](Mati & mat, std::string path, int cutoff) {
			std::ifstream ifile;
			int i, j;

			FOPEN(ifile, path);
			mat = Mati::Zero(this->m_rows + 1, this->m_bins);
			for (i = 0; i < this->m_rows + 1; i++) {
				for (j = 0; j < this->m_bins; j++) {
					ifile >> mat(i, j);
					if (j < cutoff) mat(i, j) = 0;
				}
			}
			FCLOSE(ifile);
		};

		auto bar = [&ifile](auto && mat, int L, int M, int N) {
			int i, j, l, m, n, a, b;

			for (i = 0; i < 4; i++) {
				for (j = 0; j < 4; j++) {
					ifile >> a >> b;
					for (l = 0; l < L; l++) {
						for (m = 0; m < M; m++) {
							for (n = 0; n < N; n++) {
								ifile >> mat[i][j](l, m, n);
							}
						}
					}
				}
			}
		};

		foo(m_counts_stacking, Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/stacking.counts", 0);
		foo(m_counts_pairing, Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/pairing.counts", 0);
		foo(m_counts_wc, Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/wc.counts", 0);
		foo(m_counts_nwc, Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/nwc.counts", 0);

		FOPEN(ifile, Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/pars.txt");
		bar(m_counts_bp, 40, 40, 8);
		bar(m_counts_st53, 20, 20, 20);
		bar(m_counts_st35, 20, 20, 20);
		FCLOSE(ifile);
	}

	void Score::read_pars() {
		Par par(Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/pars");
		if (par.has("bond_len_std")) {
			for (auto && i : par.getv("bond_len_std")) {
				m_bond_len_std.push_back(JN_DBL(i));
			}
		}
		if (par.has("bond_angle_std")) {
			for (auto && i : par.getv("bond_angle_std")) {
				m_bond_angle_std.push_back(JN_DBL(i));
			}
		}
		if (par.has("bond_dihedral_std")) {
			for (auto && i : par.getv("bond_dihedral_std")) {
				m_bond_dihedral_std.push_back(JN_DBL(i));
			}
		}
		par.set(m_cutoff_stacking, "cutoff_stacking");
		par.set(m_cutoff_pairing, "cutoff_pairing");
	}

	void Score::set_freqs() {
		auto foo = [](auto && counts, auto && weights, auto && freqs, int L, int M, int N) {
			int i, j, l, m, n;
			float d, sum;

			sum = 0;
			for (i = 0; i < 4; i++) {
				for (j = 0; j < 4; j++) {
					d = float(counts[i][j].sum());
					weights[i][j] = d;
					sum += d;
					for (l = 0; l < L; l++) {
						for (m = 0; m < M; m++) {
							for (n = 0; n < N; n++) {
								freqs[i][j](l, m, n) = counts[i][j](l, m, n) / d;
							}
						}
					}
				}
			}
			for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) weights[i][j] /= sum;
		};

		auto bar = [this](Mat & m, const Mati & c) {
			Vec v;
			int i, j, sum;

			m = Mat::Zero(m_rows, m_bins);
			for (i = 0; i < m_rows; i++) {
				sum = c.row(i).sum();
				for (j = 0; j < m_bins; j++) {
					m(i, j) = (sum < m_bins * 3 ? 0 : (double(c(i, j)) / sum));
				}
			}

			v = Vec::Zero(m_bins);
			sum = c.row(m_rows).sum();
			for (i = 0; i < m_bins; i++) {
				v[i] = double(c(m_rows, i)) / sum;
			}

			for (i = 0; i < m_rows; i++) {
				for (j = 0; j < m_bins; j++) {
					m(i, j) = (v[j] > 0.0003 ? (m(i, j) / v[j]) : 0);
				}
			}
		};

		bar(m_freqs_stacking, m_counts_stacking);
		bar(m_freqs_pairing, m_counts_pairing);
		bar(m_freqs_wc, m_counts_wc);
		bar(m_freqs_nwc, m_counts_nwc);

		m_counts_stacking = Mati::Zero(m_rows, m_bins);
		m_counts_pairing = Mati::Zero(m_rows, m_bins);
		m_counts_wc = Mati::Zero(m_rows, m_bins);
		m_counts_nwc = Mati::Zero(m_rows, m_bins);

		foo(m_counts_bp,   m_weights_bp,   m_freqs_bp,   40, 40, 8 );
		foo(m_counts_st53, m_weights_st53, m_freqs_st53, 20, 20, 20);
		foo(m_counts_st35, m_weights_st35, m_freqs_st35, 20, 20, 20);

	}

	double Score::en_stacking(const Residue &r1, const Residue &r2) {
		return this->en_bp(r1, r2).m_en_stacking;
	}

	double Score::en_pairing(const Residue &r1, const Residue &r2) {
		return this->en_bp(r1, r2).m_en_pairing;
	}

	bool Score::in_base(int n) {
		return true;
	}

	ScoreBase &Score::en_bp(const Residue &r1, const Residue &r2) {
		int t1, t2;
		const Residue *p1, *p2;
		Residue *temp1, *temp2;
		double d;

		auto bar = [this](auto && f, auto && w, int t1, int t2, auto && v, int l, int m, int n) -> double {
			if (std::fabs(v[0]) < l && std::fabs(v[1]) < m && std::fabs(v[2]) < n) {
				double d = f[t1][t2](int((v[0] + l) / 0.5), int((v[1] + m) / 0.5), int((v[2] + n) / 0.5));
				if (d == 0) {
					return 0.0;
				}
				else {
					double s = -std::log(d) /** w[t1][t2]*/;
					return (s > 18 ? 0 : w[t1][t2] * (s - 18));
					//return s;
				}
			}
			else {
				return 0.0;
			}
		};

		if (m_cg->is_cg(r1) && m_cg->is_cg(r2)) {
			temp1 = NULL;
			temp2 = NULL;
			p1 = &r1;
			p2 = &r2;
		}
		else {
			temp1 = new Residue(m_cg->to_cg(r1));
			temp2 = new Residue(m_cg->to_cg(r2));
			p1 = temp1;
			p2 = temp2;
		}

		t1 = pdb::res_type(r1.name);
		t2 = pdb::res_type(r2.name);
		m_en_stacking = 0;
		m_en_pairing = 0;
		m_en_vdw = 0;
		m_en_wc = 0;
		m_en_nwc = 0;
		ParBp parbp(*p1, *p2);
		//m_en_pairing += foo(t1, t2, parbp.o21_) + foo(t2, t1, parbp.o12_);
		m_en_pairing += bar(m_freqs_bp, m_weights_bp, t1, t2, parbp.o21_, 10, 10, 2);
		m_en_pairing += bar(m_freqs_bp, m_weights_bp, t2, t1, parbp.o12_, 10, 10, 2);
		//if (m_en_pairing > 1.5) LOG << ">" << m_en_pairing << "\n" << parbp.o21_ << "\n" << parbp.o12_ << "\n" << std::endl;
		m_en_stacking += bar(m_freqs_st53, m_weights_st53, t1, t2, parbp.o21_, 5, 5, 5);
		m_en_stacking += bar(m_freqs_st35, m_weights_st35, t2, t1, parbp.o12_, 5, 5, 5);

		d = geom::distance(p1->at(0), p2->at(0));
		if (d < 20 && d > 10) m_en_vdw += -25 + square(d - 15);
		//for (i = 0; i < m_res_size; i++) {
		//	if (!in_base(i)) continue;
		//	for (j = 0; j < m_res_size; j++) {
		//		if (!in_base(j)) continue;
		//		d = geom::distance(p1->at(i), p2->at(j));
		//		if (d < m_cutoff) {
		//			a = t1 * m_res_size + i;
		//			b = t2 * m_res_size + j;
		//			if (parbp.is_stacked()) {
		//				f = m_freqs_stacking(a * m_num_types + b, int(d / m_bin));
		//				m_en_stacking += (f == 0 ? 0 : -std::log(f));
		//			}
		//			if (parbp.is_wc()) {
		//				f = m_freqs_wc(a * m_num_types + b, int(d / m_bin));
		//				m_en_wc += (f == 0 ? 0 : -std::log(f));
		//			}
		//			if (parbp.is_nwc()) {
		//				f = m_freqs_nwc(a * m_num_types + b, int(d / m_bin));
		//				m_en_nwc += (f == 0 ? 0 : -std::log(f));
		//			}
		//			f = m_freqs_pairing(a * m_num_types + b, int(d / m_bin));
		//			m_en_pairing += (f == 0 ? 0 : -std::log(f));
		//		}
		//	}
		//}
		if (temp1 != NULL) delete temp1;
		if (temp2 != NULL) delete temp2;
		//m_en_stacking = (m_en_stacking > m_cutoff_stacking ? 0 : m_en_stacking);
		//m_en_stacking = (m_en_stacking > );
		//if (m_en_stacking > 0) m_en_stacking -= 3;
		//if (m_en_pairing > 0) m_en_pairing -= 5;
		return *this;
	}

	void Score::print_counts(std::ostream & stream) const {
		auto foo = [&stream](auto && m) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					stream << i << ' ' << j << std::endl;
					stream << m[i][j] << std::endl;
				}
			}
		};

		foo(m_counts_bp);
		foo(m_counts_st53);
		foo(m_counts_st35);

		stream << "Stacking counts:\n" << m_counts_stacking << std::endl;
		stream << "Pairing counts:\n" << m_counts_pairing << std::endl;
		stream << "WC counts:\n" << m_counts_wc << std::endl;
		stream << "nWC counts:\n" << m_counts_nwc << std::endl;
	}

	void Score::print_freqs(std::ostream & stream) const {
		auto foo = [&stream](auto && m) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					stream << i << ' ' << j << std::endl;
					stream << m[i][j] << std::endl;
				}
			}
		};

		foo(m_freqs_bp);
		foo(m_freqs_st53);
		foo(m_freqs_st35);

		stream << "Stacking freqs:\n" << m_freqs_stacking << std::endl;
		stream << "Pairing freqs:\n" << m_freqs_pairing << std::endl;
	}

	double Score::en_len(const Chain &c, int beg) {
		double d;

		d = geom::distance(c[beg][0], c[beg + 1][0]);
		return square(d - m_bond_len_std[0]);
	}

	double Score::en_ang(const Chain &c, int beg) {
		double d;

		d = geom::angle(c[beg][0], c[beg + 1][0], c[beg + 2][0]);
		return square(d - m_bond_angle_std[0]);
	}

	double Score::en_dih(const Chain &c, int beg) {
		double d;

		d = geom::dihedral(c[beg][0], c[beg + 1][0], c[beg + 2][0], c[beg + 3][0]);
		d = d - m_bond_dihedral_std[0];
		d = 3.3 - 4 * std::cos(d) + std::cos(2 * d) - 0.44 * std::cos(3 * d);
		return d;
	}

	double Score::en_crash(const Residue &r1, const Residue &r2) {
		int i, j;
		double d, e;

		e = 0;
		for (i = 0; i < m_res_size; i++) {
			for (j = 0; j < m_res_size; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (i == 0 && j == 0) {
					if (d < 5) {
						e += square(d - 5);
					}
				}
				else if (i == 1 && j == 1) {
					if (d < 5) {
						e += square(d - 5);
					}
				}
				else if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
					if (d < 6.5) {
						e += square(d - 6.5);
					}
				}
				else if (d < 3.5) {
					e += square(d - 3.5);
				}
			}
		}
		return e;
	}

} // namespace jian

