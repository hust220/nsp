#include "ScorePsb.hpp"

namespace jian {

	REG_SCORER("psb", Score<CGpsb>);

	void Score<CGpsb>::init() {
		m_bins = int(std::ceil(m_cutoff / m_bin));
		read_counts();
		set_freqs();
	}

	//void Score<CGpsb>::train(const std::string &s) {
	//	m_counts_stacking = Mati::Zero(144, m_bins);
	//	m_counts_pairing = Mati::Zero(144, m_bins);
	//	EACH_SPLIT_LINE(s.c_str(), " ",
	//		std::cout << F[0] << std::endl;
	//		Chain chain;
	//		chain_read_model(chain, F[0]);
	//		update_counts(CGpsb::chain(chain));
	//	);
	//	std::cout << m_counts_stacking << std::endl;
	//	std::cout << m_counts_pairing << std::endl;
	//}

	void Score<CGpsb>::run(const Chain & c) {
		int i, j, l;
		Chain chain;

		chain = CGpsb::chain(c);
		l = chain.size();
		m_score = 0;
		for (i = 0; i < l; i++) {
			m_score += en_stacking(chain[i], chain[j]);
			for (j = i + 2; j < l; j++) {
				m_score += en_pairing(chain[i], chain[j]);
			}
		}
	}

	void Score<CGpsb>::train(const Chain &c) {
		unsigned i, j, len;
		Chain chain;

		chain = CGpsb::chain(c);
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

	void Score<CGpsb>::set_indices() {
		unsigned i, len;

		len = m_chain->size();
		m_indices.resize(len);
		for (i = 0; i < len; i++) {
			m_indices[i] = m_map[m_chain->at(i).name];
		}
	}

	void Score<CGpsb>::update_counts_stacking(int n1, int n2) {
		unsigned i, j;
		double d;
		int n, a = m_indices[n1], b = m_indices[n2];
		const Residue &r1 = m_chain->at(n1);
		const Residue &r2 = m_chain->at(n2);
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (d < 20) {
					n = (a * 3 + i) * 12 + (b * 3 + j);
					m_counts_stacking(n, int(d / m_bin))++;
				}
			}
		}
	}

	void Score<CGpsb>::update_counts_pairing(int n1, int n2) {
		unsigned i, j;
		double d;
		int n, a = m_indices[n1], b = m_indices[n2];
		const Residue &r1 = m_chain->at(n1);
		const Residue &r2 = m_chain->at(n2);
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (d < 20) {
					n = (a * 3 + i) * 12 + (b * 3 + j);
					m_counts_pairing(n, int(d / m_bin))++;
				}
			}
		}
	}

	void Score<CGpsb>::read_counts() {
		m_counts_stacking = Mati::Zero(144, m_bins);
		m_counts_pairing = Mati::Zero(144, m_bins);
		std::string name = Env::lib() + "/RNA/pars/scoring/score_psb/counts";
		std::ifstream ifile(name.c_str());
		for (int i = 0; i < 144; i++) {
			for (int j = 0; j < m_bins; j++) {
				ifile >> m_counts_stacking(i, j);
			}
		}
		for (int i = 0; i < 144; i++) {
			for (int j = 0; j < m_bins; j++) {
				ifile >> m_counts_pairing(i, j);
			}
		}
		ifile.close();
	}

	void Score<CGpsb>::set_freqs() {
		m_freqs_stacking = Mat::Zero(144, m_bins);
		m_freqs_pairing = Mat::Zero(144, m_bins);
		Vec vs = Vec::Zero(m_bins);
		Vec vp = Vec::Zero(m_bins);
		int n, sum_s, sum_p, total_sum_s = 0, total_sum_p = 0;
		for (int i = 0; i < 144; i++) {
			sum_s = 0;
			sum_p = 0;
			for (int j = 0; j < m_bins; j++) {
				n = m_counts_stacking(i, j);
				sum_s += n;
				total_sum_s += n;
				vs[j] += n;
				n = m_counts_pairing(i, j);
				sum_p += n;
				total_sum_p += n;
				vp[j] += n;
			}
			for (int j = 0; j < m_bins; j++) {
				m_freqs_stacking(i, j) = m_counts_stacking(i, j) / double(sum_s);
				m_freqs_pairing(i, j) = m_counts_pairing(i, j) / double(sum_p);
			}
		}
		for (int i = 0; i < m_bins; i++) {
			vs[i] /= double(total_sum_s);
			vp[i] /= double(total_sum_p);
		}
		for (int i = 0; i < 144; i++) {
			for (int j = 0; j < m_bins; j++) {
				if (vs[j] > 0.0003) {
					m_freqs_stacking(i, j) /= vs[j];
				}
				else {
					m_freqs_stacking(i, j) = 0;
				}
				if (vp[j] > 0.0003) {
					m_freqs_pairing(i, j) /= vp[j];
				}
				else {
					m_freqs_pairing(i, j) = 0;
				}
			}
		}
		//        std::cout << m_freqs_stacking << std::endl;
	}

	double Score<CGpsb>::en_stacking(const Residue &r1, const Residue &r2) {
		double e = 0, d, f;
		int a, b, t1 = m_map[r1.name], t2 = m_map[r2.name];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (d < 20) {
					a = t1 * 3 + i;
					b = t2 * 3 + j;
					f = m_freqs_stacking(a * 12 + b, int(d / m_bin));
					if (f != 0) {
						e += -std::log(f);
					}
					else {
						e += 0;
					}
				}
			}
		}
		return e;
	}

	double Score<CGpsb>::en_pairing(const Residue &r1, const Residue &r2) {
		double e = 0, d, f;
		int a, b, t1 = m_map[r1.name], t2 = m_map[r2.name];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				d = geom::distance(r1[i], r2[j]);
				if (d < 20) {
					a = t1 * 3 + i;
					b = t2 * 3 + j;
					f = m_freqs_pairing(a * 12 + b, int(d / m_bin));
					if (f != 0) {
						e += -std::log(f);
					}
					else {
						e += 0;
					}
				}
			}
		}
		return e;
	}

	void Score<CGpsb>::print_freqs_stacking(std::ostream & stream) const {
		stream << m_freqs_stacking << std::endl;
	}

	void Score<CGpsb>::print_freqs_pairing(std::ostream & stream) const {
		stream << m_freqs_pairing << std::endl;
	}

	void Score<CGpsb>::print_counts_stacking(std::ostream & stream) const {
		stream << m_counts_stacking << std::endl;
	}

	void Score<CGpsb>::print_counts_pairing(std::ostream & stream) const {
		stream << m_counts_pairing << std::endl;
	}

	void Score<CGpsb>::print_counts(std::ostream & stream) const {
		stream << "Stacking counts: " << std::endl;
		print_counts_stacking(stream);
		stream << "Pairing counts: " << std::endl;
		print_counts_pairing(stream);
	}

} // namespace jian

