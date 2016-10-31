#include "Score.hpp"
#include "../utils/Par.hpp"

namespace jian {

	namespace cg1p { REG_SCORER("1p", Score); }
	namespace cgpsb { REG_SCORER("psb", Score); }
	namespace cg6p { REG_SCORER("6p", Score); }

	Score::Score(std::string cg) {
		m_cg = CG::fac_t::create(cg);
	}

	void Score::init() {
		m_bin = 0.3;
		m_cutoff = 20;
		m_bins = int(std::ceil(m_cutoff / m_bin));
		m_res_size = m_cg->res_size();
		m_num_types = m_res_size * 4;
		m_rows = m_num_types * m_num_types;
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

		chain = m_cg->to_cg(c);
		m_chain = &chain;
		set_indices();
		len = chain.size();
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

	void Score::update_counts_pairing(int n1, int n2) {
		int i, j;
		double d;
		int n, a = m_indices[n1], b = m_indices[n2];
		const Residue &r1 = m_chain->at(n1);
		const Residue &r2 = m_chain->at(n2);
		for (i = 0; i < m_res_size; i++) {
			for (j = 0; j < m_res_size; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (d < 20) {
					n = (a * m_res_size + i) * m_num_types + (b * m_res_size + j);
					m_counts_pairing(n, int(d / m_bin))++;
				}
			}
		}
	}

	void Score::read_counts(Mati & m, std::string path, int cutoff) {
		std::ifstream ifile;
		int i, j;

		FOPEN(ifile, path);
		m = Mati::Zero(m_rows + 1, m_bins);
		for (i = 0; i < m_rows + 1; i++) {
			for (j = 0; j < m_bins; j++) {
				ifile >> m(i, j);
				if (j < cutoff) m(i, j) = 0;
			}
		}
		FCLOSE(ifile);
	}

	void Score::read_counts() {
		read_counts(m_counts_stacking, Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/stacking.counts", 0);
		read_counts(m_counts_pairing, Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/pairing.counts", 0);
	}

	void Score::read_pars() {
		Par par(Env::lib() + "/RNA/pars/scoring/score_" + m_cg->m_cg + "/pars");
		par.set(m_bond_len_std, "bond_len_std");
		par.set(m_bond_angle_std, "bond_angle_std");
		par.set(m_bond_dihedral_std, "bond_dihedral_std");
	}

	void Score::set_freqs(Mat & m, const Mati & c) {
		Vec v;
		int i, j, sum;

		m = Mat::Zero(m_rows, m_bins);
		for (i = 0; i < m_rows; i++) {
			sum = c.row(i).sum();
			for (j = 0; j < m_bins; j++) {
				m(i, j) = (sum == 0 ? 0 : (double(c(i, j)) / sum));
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
	}

	void Score::set_freqs() {
		int i, j;

		for (i = 0; i < m_rows; i++) for (j = 0; j < m_bins; j++) m_counts_pairing(i, j) += m_counts_stacking(i, j);
		set_freqs(m_freqs_stacking, m_counts_stacking);
		set_freqs(m_freqs_pairing, m_counts_pairing);
		m_counts_stacking = Mati::Zero(m_rows, m_bins);
		m_counts_pairing = Mati::Zero(m_rows, m_bins);

	}

	double Score::en_stacking(const Residue &r1, const Residue &r2) {
		double e, d, f;
		int a, b, t1, t2;
		const Residue *p1, *p2;
		Residue *temp1, *temp2;

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

		e = 0;
		t1 = m_map[p1->name];
		t2 = m_map[p1->name];
		for (int i = 0; i < m_res_size; i++) {
			for (int j = 0; j < m_res_size; j++) {
				d = geom::distance(p1->at(i), p2->at(j));
				if (d < 20) {
					a = t1 * m_res_size + i;
					b = t2 * m_res_size + j;
					f = m_freqs_stacking(a * m_num_types + b, int(d / m_bin));
					if (f != 0) {
						e += -std::log(f);
					}
					else {
						e += 0;
					}
				}
			}
		}
		if (temp1 != NULL) delete temp1;
		if (temp2 != NULL) delete temp2;
		return (e > 0 ? 0 : e);
	}

	double Score::en_pairing(const Residue &r1, const Residue &r2) {
		double e, d, f;
		int i, j, a, b, t1, t2;
		const Residue *p1, *p2;
		Residue *temp1, *temp2;

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

		t1 = m_map[r1.name];
		t2 = m_map[r2.name];
		e = 0;
		for (i = 0; i < m_res_size; i++) {
			for (j = 0; j < m_res_size; j++) {
				d = geom::distance(p1->at(i), p2->at(j));
				if (d < m_cutoff) {
					a = t1 * m_res_size + i;
					b = t2 * m_res_size + j;
					f = m_freqs_pairing(a * m_num_types + b, int(d / m_bin));
					e += (f == 0 ? 0 : -std::log(f));
				}
			}
		}
		if (temp1 != NULL) delete temp1;
		if (temp2 != NULL) delete temp2;
		return (e > 0 ? 0 : e);
	}

	void Score::print_counts(std::ostream & stream, const Mati & c) const {
		stream << c << std::endl;
	}

	void Score::print_counts(std::ostream & stream) const {
		stream << "Stacking counts: " << std::endl;
		print_counts(stream, m_counts_stacking);
		stream << "Pairing counts: " << std::endl;
		print_counts(stream, m_counts_pairing);
	}

	void Score::print_freqs(std::ostream & stream, const Mat & f) const {
		stream << f << std::endl;
	}

	void Score::print_freqs(std::ostream & stream) const {
		stream << "Stacking freqs: " << std::endl;
		print_freqs(stream, m_freqs_stacking);
		stream << "Pairing freqs: " << std::endl;
		print_freqs(stream, m_freqs_pairing);
	}

	double Score::en_len(const Residue &r1, const Residue &r2) {
		double d;

		if (m_cg->m_cg == "psb") {
			d = geom::distance(r1[0], r2[0]);
			return square(d - m_bond_len_std);
		}
		else {
			d = geom::distance(r1[2], r2[0]);
			return square(d - 3.1);
		}
	}

	double Score::en_ang(const Residue &r1, const Residue &r2, const Residue &r3) {
		double d;

		d = geom::angle(r1[0], r2[0], r3[0]);
		return square(d - m_bond_angle_std);
	}

	double Score::en_dih(const Residue &r1, const Residue &r2, const Residue &r3, const Residue &r4) {
		double d;

		d = geom::dihedral(r1[0], r2[0], r3[0], r4[0]);
		d = d - m_bond_dihedral_std;
		d = 3.3 - 4 * std::cos(d) + std::cos(2 * d) - 0.44 * std::cos(3 * d);
		return d;
	}

	double Score::en_crash_6p(const Residue &r1, const Residue &r2) {
		int i, j;
		double d, e;

		e = 0;
		for (i = 0; i < m_res_size; i++) {
			for (j = 0; j < m_res_size; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (d < 3.5) {
					e += square(d - 3.5);
				}
			}
		}
		return e;
	}

	double Score::en_crash_psb(const Residue &r1, const Residue &r2) {
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

	double Score::en_crash(const Residue &r1, const Residue &r2) {
		if (m_cg->m_cg == "psb") {
			return en_crash_psb(r1, r2);
		}
		else if (m_cg->m_cg == "6p") {
			return en_crash_6p(r1, r2);
		}
		else {
			throw "jian::Score::en_crash error!";
		}
	}

} // namespace jian

