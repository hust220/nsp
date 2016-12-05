#include "../pdb.hpp"
#include "../geom.hpp"
#include "DistAnal.hpp"

BEGIN_JN

#define FREE_COORDS do {\
	for (auto && p : m_coords) {\
		delete [] p;\
	}\
} while (0)

#define PRINT_MAT3(a, stream) do { \
	int i, j, k; \
	for (i = 0; i < 85; i++) {\
		for (j = 0; j < 85; j++) {\
			for (k = 0; k < bins; k++) {\
				stream << a[(i * 85 + j) * bins + k] << ' ';\
			}\
			stream << std::endl;\
		}\
	}\
} while (0)


	DistAnal & DistAnal::init(double interval, double cutoff) {
		this->interval = interval;
		this->cutoff = cutoff;
		bins = int(ceil(cutoff / interval));
		m_freqs.resize(85 * 85 * bins);
		m_counts.resize(85 * 85 * bins);
		return *this;
	}

	DistAnal::~DistAnal() {
		FREE_COORDS;
	}

	inline int DistAnal::res_type(S name) {
		return (name == "A" ? 1 : (name == "U" ? 2 : (name == "G" ? 3 : (name == "C" ? 4 : (name == "T" ? 2 : -1)))));
	}

	inline int DistAnal::atom_type(const Residue &r, int k) {
		int t = res_type(r.name);
		return (r[0].name == "P" ? 0 : 3) + (t == 1 ? k : (t == 2 ? 22 + k : (t == 3 ? 42 + k : 65 + k)));
	}

	void DistAnal::read_mol(const Chain &chain) {
		len = chain.size();
		num.resize(len);
		type.resize(len);
		size_nt.resize(len);
		FREE_COORDS;
		m_coords.resize(len);

		Point p;
		int m = 0, n = 0;
		for (auto && res : chain) {
			m++;
			size_nt[n] = res.size();
			type[n] = res_type(res.name);
			m_coords[n] = new Point[size_nt[n]];
			int k = 0;
			for (auto && atom : res) {
				for (int l = 0; l < 3; l++) m_coords[n][k][l] = atom[l];
				m_coords[n][k][3] = atom_type(res, k);
				//m_coords[n][k][3] = (res[0].name == "P" ? 0 : 3);
				//m_coords[n][k][3] += (type[n] == 1 ? k : (type[n] == 2 ? 22 + k : (type[n] == 3 ? 42 + k : 65 + k)));
				if (atom.name == "O5*" && n != 0) {
					double dist = geom::distance(m_coords[n][k], p);
					if (dist > 4) m += 10000;
				}
				else if (atom.name == "O3*") {
					for (int l = 0; l < 3; l++) p[l] = m_coords[n][k][l];
				}
				k++;
			}
			num[n] = m;
			n++;
		}
	}

	DistAnal & DistAnal::train(const Chain & c) {
		int i, j, k, l, type1, type2;
		double temp;

		read_mol(c);
		for (i = 0; i < len; i++) {
			for (j = i + 1; j < len; j++) {
				for (k = 0; k < size_nt[i]; k++) {
					for (l = 0; l < size_nt[j]; l++) {
						type1 = static_cast<int>(m_coords[i][k][3]);
						type2 = static_cast<int>(m_coords[j][l][3]);
						if (num[j] - num[i] == 1 && !in_base(type1) && !in_base(type2)) continue;
						temp = geom::distance(m_coords[i][k], m_coords[j][l]);
						if (temp >= cutoff) continue;
						m_counts[static_cast<int>(m_coords[i][k][3] * 85 + m_coords[j][l][3]) * bins +
							     static_cast<int>(temp / interval)]++;
						m_counts[static_cast<int>(m_coords[j][l][3] * 85 + m_coords[i][k][3]) * bins +
							     static_cast<int>(temp / interval)]++;
					}
				}
			}
		}
		return *this;
	}

	bool DistAnal::in_base(int type) {
		return (type > 11 && type < 22) ||
			   (type > 33 && type < 42) ||
			   (type > 54 && type < 62) ||
			   (type > 74 && type < 85);
	}

	DistAnal & DistAnal::run(const Chain &chain) {
		int i, j, k, l, type1, type2;
		double a, b, temp;

		read_mol(chain);
		score = 0;
		for (i = 0; i < len; i++) {
			for (j = i + 1; j < len; j++) {
				for (k = 0; k < size_nt[i]; k++) {
					for (l = 0; l < size_nt[j]; l++) {
						type1 = int(m_coords[i][k][3]);
						type2 = int(m_coords[j][l][3]);
						if (num[j] - num[i] == 1 && !in_base(type1) && !in_base(type2)) continue;
						temp = geom::distance(m_coords[i][k], m_coords[j][l]);
						if (temp >= cutoff) continue;
						a = m_freqs[(type1 * 85 + type2) * bins + int(temp / interval)];
						b = m_freqs[(type2 * 85 + type1) * bins + int(temp / interval)];
						score += (a == 0 ? penalty : (-log(a) * ((in_base(type1) && in_base(type2)) ? 2.5 : 1)));
						score += (b == 0 ? penalty : (-log(b) * ((in_base(type1) && in_base(type2)) ? 2.5 : 1)));
					}
				}
			}
		}
		score = score / (len * (len - 1));
		return *this;
	}

	void DistAnal::read_freqs(S filename) {
		int i, j, k;
		std::ifstream ifile;

		FOPEN(ifile, filename);
		for (i = 0; i < 85; i++) {
			for (j = 0; j < 85; j++) {
				for (k = 0; k < bins; k++) {
					ifile >> m_freqs[(i * 85 + j) * bins + k];
				}
			}
		}
		FCLOSE(ifile);
	}

	void DistAnal::print_freqs(std::ostream & stream) const {
		PRINT_MAT3(m_freqs, stream);
	}

	void DistAnal::print_counts(std::ostream & stream) const {
		PRINT_MAT3(m_counts, stream);
	}

	double DistAnal::en_stacking(const Residue &r1, const Residue &r2) {
		int n1, n2, t1, t2;
		double d, a, b, e;

		e = 0;
		n1 = 0;
		for (auto && a1 : r1) {
			t1 = atom_type(r1, n1);
			n2 = 0;
			for (auto && a2 : r2) {
				t2 = atom_type(r2, n2);
				if (!in_base(t1) && !in_base(t2)) continue;
				d = geom::distance(a1, a2);
				if (d >= cutoff) continue;
				a = m_freqs[(t1 * 85 + t2) * bins + int(d / interval)];
				b = m_freqs[(t2 * 85 + t1) * bins + int(d / interval)];
				e += (a == 0 ? penalty : (-log(a) * ((in_base(t1) && in_base(t2)) ? 2.5 : 1)));
				e += (b == 0 ? penalty : (-log(b) * ((in_base(t1) && in_base(t2)) ? 2.5 : 1)));
				n2++;
			}
			n1++;
		}
		return e;
	}

	double DistAnal::en_pairing(const Residue &r1, const Residue &r2) {
		int n1, n2, t1, t2;
		double d, a, b, e;

		e = 0;
		n1 = 0;
		for (auto && a1 : r1) {
			t1 = atom_type(r1, n1);
			n2 = 0;
			for (auto && a2 : r2) {
				t2 = atom_type(r2, n2);
				d = geom::distance(a1, a2);
				if (d >= cutoff) continue;
				a = m_freqs[(t1 * 85 + t2) * bins + int(d / interval)];
				b = m_freqs[(t2 * 85 + t1) * bins + int(d / interval)];
				e += (a == 0 ? penalty : (-log(a) * ((in_base(t1) && in_base(t2)) ? 2.5 : 1)));
				e += (b == 0 ? penalty : (-log(b) * ((in_base(t1) && in_base(t2)) ? 2.5 : 1)));
				n2++;
			}
			n1++;
		}
		return e;
	}


END_JN

