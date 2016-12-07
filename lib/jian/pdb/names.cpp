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

			for (auto &&it : FileLines(filename, "#: []")) {
				int l = size(it.arr);
				if (l > 0) {
					//std::cout << F[0] << std::endl;
					if (it.arr[0] == "mol_type") {
						//map_names[F[1]] = Names{};
						p = map_names[it.arr[1]];
					}
					else if (it.arr[0] == "phos") {
						for (i = 1; i < l; i++) {
							p->atoms_phos.push_back(it.arr[i]);
							p->atoms_bb.push_back(it.arr[i]);
						}
					}
					else if (it.arr[0] == "sugar") {
						for (i = 1; i < l; i++) {
							p->atoms_sugar.push_back(it.arr[i]);
							p->atoms_bb.push_back(it.arr[i]);
						}
					}
					else if (it.arr[0] == "bb") {
						for (i = 1; i < l; i++) {
							p->atoms_bb.push_back(it.arr[i]);
						}
					}
					else {
						jian::tokenize(it.line, v, ":");
						//for (auto && i : v) std::cout << i << ' '; std::cout << std::endl;
						jian::tokenize(std::string(v[0]), v, " []");
						for (i = 1; i < v.size(); i++) p->alias[it.arr[0]].push_back(v[i]);

						p->res.push_back(it.arr[0]);
						for (auto && s : p->atoms_phos) p->atoms_res[it.arr[0]].push_back(s);
						for (auto && s : p->atoms_sugar) p->atoms_res[it.arr[0]].push_back(s);
						for (i = v.size(); i < it.arr.size(); i++) {
							p->atoms_base[it.arr[0]].push_back(it.arr[i]);
							p->atoms_res[it.arr[0]].push_back(it.arr[i]);
						}
					}
				}
				else {
					// pass
				}
			}
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

