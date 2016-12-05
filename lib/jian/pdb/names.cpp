#include <memory>
#include <algorithm>
#include <regex>
#include "../utils/file.hpp"
#include "../utils/Env.hpp"
#include "names.hpp"

BEGIN_JN
	namespace pdb {
		Names::Names() {}

		void Names::print_names(const Names & names) {
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

		void Names::init_map(Names::map_Names &map_names) {
			std::shared_ptr<Names> p;
			S filename = Env::lib() + "/RNA/pars/pdb/names";
			int i;
			jian::tokenize_v v;

			BEGIN_READ_FILE(filename, "#: []") {
				if (F.size() > 0) {
					//std::cout << F[0] << std::endl;
					if (F[0] == "mol_type") {
						//map_names[F[1]] = Names{};
						p = map_names[F[1]];
					}
					else if (F[0] == "phos") {
						for (i = 1; i < F.size(); i++) {
							p->atoms_phos.push_back(F[i]);
							p->atoms_bb.push_back(F[i]);
						}
					}
					else if (F[0] == "sugar") {
						for (i = 1; i < F.size(); i++) {
							p->atoms_sugar.push_back(F[i]);
							p->atoms_bb.push_back(F[i]);
						}
					}
					else if (F[0] == "bb") {
						for (i = 1; i < F.size(); i++) {
							p->atoms_bb.push_back(F[i]);
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
		}

		const Names & Names::instance(S mol_type) {
			static std::map<std::string, std::shared_ptr<Names>> map;
			if (map.empty()) {
				map["RNA"].reset(new Names);
				map["DNA"].reset(new Names);
				map["protein"].reset(new Names);
				init_map(map);
			}
			if (map.find(mol_type) == map.end()) {
				return *map["RNA"];
			}
			else {
				return *map[mol_type];
			}
		}

	} // namespace pdb
END_JN

