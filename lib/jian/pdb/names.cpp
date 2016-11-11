#include <algorithm>
#include <regex>
#include "../utils/file.hpp"
#include "../utils/Env.hpp"
#include "names.hpp"

namespace jian {
	namespace pdb {
		using map_Names = std::map<std::string, Names>;

		namespace names_detail {
			void print_names(const Names & names) {
				for (auto && i : names.alias) {
					std::cout << i.first << ' ';
					for (auto && j : i.second) std::cout << j << ' '; std::cout << std::endl;
				}
				for (auto && i : names.atoms_base) {
					std::cout << i.first << ' ';
					for (auto && j : i.second) std::cout << j << ' '; std::cout << std::endl;
				}
				for (auto && i : names.atoms_res) {
					std::cout << i.first << ' ';
					for (auto && j : i.second) std::cout << j << ' '; std::cout << std::endl;
				}
			}

			map_Names init_names() {
				map_Names map_names;
				Names *p;
				std::string filename = Env::lib() + "/RNA/pars/pdb/names";
				int i;
				jian::tokenize_v v;

				BEGIN_READ_FILE(filename, "#: []") {
					if (F.size() > 0) {
						//std::cout << F[0] << std::endl;
						if (F[0] == "mol_type") {
							map_names[F[1]] = Names{};
							//map_names.emplace(F[1]);
							p = &(map_names[F[1]]);
						}
						else if (F[0] == "phos") {
							for (i = 1; i < F.size(); i++) {
								p->atoms_phos.push_back(F[i]);
							}
						}
						else if (F[0] == "sugar") {
							for (i = 1; i < F.size(); i++) {
								p->atoms_sugar.push_back(F[i]);
							}
						}
						else {
							jian::tokenize(L, v, ":");
							//for (auto && i : v) std::cout << i << ' '; std::cout << std::endl;
							jian::tokenize(std::string(v[0]), v, " []");
							for (i = 1; i < v.size(); i++) p->alias[F[0]].push_back(v[i]);

							p->res.push_back(F[0]);
							for (auto && s : p->atoms_phos) p->atoms_res[F[0]].push_back(s);
							for (auto && s : p->atoms_sugar) p->atoms_res[F[0]].push_back(s);
							for (i = v.size(); i < F.size(); i++) {
								p->atoms_base[F[0]].push_back(F[i]);
								p->atoms_res[F[0]].push_back(F[i]);
							}
						}
					}
					else {
						// pass
					}
				} END_READ_FILE;
				//for (auto && names : map_names) {
				//	print_names(names.second);
				//}
				return map_names;
			}

		}

		const Names & Names::instance(std::string mol_type) {
			static map_Names names = names_detail::init_names();
			if (names.find(mol_type) == names.end()) {
				return names["RNA"];
			}
			else {
				return names[mol_type];
			}
		}

	} // namespace pdb
} // namespace jian

