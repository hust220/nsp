#include "EnPsb.hpp"
#include "jian/utils/file.hpp"
#include "jian/utils/Par.hpp"

namespace jian {

	void EnPsb::init() {
		read_parameters();
		read_stacking_pairing_parameters();
	}

#define SET_MEM_EN(a) temp_par.set(PP_CAT(m_, a), PP_STRING3(PP_CAT(mc_, a)));
	void EnPsb::read_parameters() {
		Par temp_par(Env::lib() + "/RNA/pars/nuc3d/mc/mcpsb.par");
		JN_MAP(SET_MEM_EN, MEM_EN);
	}

	void EnPsb::read_stacking_pairing_parameters() {
		int i, j;
		std::string file_name;

		m_stacking_pars.resize(8, 8);
		file_name = Env::lib() + "/RNA/pars/nuc3d/mc/mcpsb.stacking.par";
		//std::cout << file_name << std::endl;
		std::ifstream ifile;
		FOPEN(ifile, file_name);
		for (i = 0; i < 8; i++) for (j = 0; j < 8; j++) ifile >> m_stacking_pars(i, j);
		ifile.close();
		m_pairing_pars.resize(3, 3);
		file_name = Env::lib() + "/RNA/pars/nuc3d/mc/mcpsb.pairing.par";
		ifile.open(file_name.c_str());
		for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) ifile >> m_pairing_pars(i, j);
		FCLOSE(ifile);
	}

#define PRINT_PARS(a) std::cout << #a << ": " << m_##a << std::endl;
	void EnPsb::print_parameters() const {
		std::cout << "parameters:" << std::endl;
		JN_MAP(PRINT_PARS, MEM_EN);
		std::cout << "stacking parameters:" << std::endl;
		std::cout << m_stacking_pars << std::endl;
		std::cout << "pairing parameters:" << std::endl;
		std::cout << m_pairing_pars << std::endl;
	}


	double EnPsb::en_stacking(const Mat &m, int a, int b) const {
		double en = 0, d;
		int x, y;
		int i, j;

		x = a * 2;
		y = b * 2;
		en = 0;
		if (m(0, 0) > 4 && m(0, 0) < 7 && m(1, 1) > 3.5 && m(1, 1) < 5.5) {
			for (i = 0; i < 2; i++) {
				for (j = 0; j < 2; j++) {
					d = square(m(i, j) - m_stacking_pars(x + i, y + j));
					if (d > 1) return 0;
					en += d;
				}
			}
			return (-40 + en) * m_stacking_weight;
		}
		else {
			return 0;
		}
	}

	double EnPsb::en_pairing(const Mat & m, int a, int b) const {
		double d, en;
		int x, y;
		unsigned i, j;

		en = 0;
		x = a + b;
		y = a*b;
		if (m(0, 0) > 14.5 && m(0, 0) < 17.5 && m(1, 1) > 9.5 && m(1, 1) < 14.5) {
			if (a == 0 || a == 2) {
				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						d = square(m(i, j) - m_pairing_pars(i, j));
						if (d > 1) return 0;
						en += d;
					}
				}
			}
			else {
				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						d = square(m(i, j) - m_pairing_pars(j, i));
						if (d > 1) return 0;
						en += d;
					}
				}
			}
			if (x == 1) {
				return (-100 + en) * 2 * m_pairing_weight;
			}
			else if (x == 5) {
				return (-100 + en) * 3 * m_pairing_weight;
			}
			else if (y == 2) {
				return (-100 + en) * 1 * m_pairing_weight;
			}
			else {
				return (-100 + en) * 0.5 * m_pairing_weight;
			}
		}
		else {
			return 0;
		}
	}

	double EnPsb::en_crash(const Residue & a, const Residue & b, Mat & arr) const {
		double d, e;
		int i, j;

		e = 0;
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				d = geom::distance(a[i], b[j]);
				arr(i, j) = d;
				if (i == 0 || j == 0) {
					if (d < 7) {
						e += square(d - 7);
					}
				}
				else if ((i == 1 || j == 1) && d < 5) {
					e += square(d - 5);
				}
				else if (d < 3.5) {
					e += square(d - 3.5);
				}
			}
		}
		return e * m_crash_weight;
	}

} // namespace jian
