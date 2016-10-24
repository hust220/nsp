#include "nsp.hpp"
#include <jian/pdb.hpp>
#include <jian/mcsm/EnPsb.hpp>

namespace jian {
	namespace {
		void set_indices(const Chain & c, std::vector<int> & v) {
			static std::vector<std::string> names{ "A", "U", "G", "C" };
			int l = c.size();
			v.resize(l);
			for (int i = 0; i < l; i++) {
				v[i] = std::distance(names.begin(), std::find(names.begin(), names.end(), c[i].name));
			}
		}

		REGISTER_NSP_COMPONENT(en_psb) {
			Chain c;
			int i, j, l;
			std::string file_name = par.get("s", "pdb", "cif");
			std::vector<int> v;
			EnPsb en_psb;
			Mat arr(3, 3);
			//double en_crash, en_stacking, en_pairing;

			chain_read_model(c, file_name);
			c = CGpsb::chain(c);
			//std::cout << c << std::endl;
			l = c.size();
			set_indices(c, v);
			for (auto && i : v) std::cout << i << ' '; std::cout << std::endl;
			en_psb.init();
			en_psb.print_parameters();

			for (i = 0; i < l; i++) {
				for (j = i + 1; j < l; j++) {
					std::cout << i << ' ' << j << ' ';
					std::cout << en_psb.en_crash(c[i], c[j], arr) << "(crash)" << std::endl;
					std::cout << arr << std::endl;
					std::cout << en_psb.en_stacking(arr, v[i], v[j]) << "(stacking) ";
					std::cout << en_psb.en_pairing(arr, v[i], v[j]) << "(pairing)" << std::endl;
				}
			}
		}
	}
} // namespace jian

