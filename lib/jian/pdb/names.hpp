#pragma once

#include <memory>
#include <string>
#include <vector>
#include <map>

BEGIN_JN
	namespace pdb {

		using names_t = std::vector<std::string>;
		using map_names_t = std::map<std::string, names_t>;

		class Names {
		public:
			using map_Names = std::map<std::string, std::shared_ptr<Names>>;

			static const Names & instance(S mol_type);

			names_t res;
			map_names_t alias;
			map_names_t atoms_base;
			map_names_t atoms_res;
			names_t atoms_phos;
			names_t atoms_sugar;
			names_t atoms_bb;

		private:
			Names();
			static void init_map(map_Names &map);
			void print_names(const Names & names);
		};

		inline int res_type(const S &res_name) {
			for (auto && mol_type : { "RNA", "DNA", "protein" }) {
				const Names &names = Names::instance(mol_type);
				auto it = std::find_if(names.res.begin(), names.res.end(), [&res_name, &names](const S &s) {
					const names_t &v = names.alias.at(s);
					return s == res_name || std::find(v.begin(), v.end(), res_name) != v.end();
				});
				if (it != names.res.end()) {
					return std::distance(names.res.begin(), it);
				}
			}
			throw std::string("unknown residue: ") + res_name;
		}


		template<typename T>
		inline T res_type(const S &res_name) {
			for (auto && mol_type : { "RNA", "DNA", "protein" }) {
				const Names &names = Names::instance(mol_type);
				auto it = std::find_if(names.res.begin(), names.res.end(), [&res_name, &names](const S &s) {
					const names_t &v = names.alias.at(s);
					return s == res_name || std::find(v.begin(), v.end(), res_name) != v.end();
				});
				if (it != names.res.end()) {
					return *it;
				}
			}
			throw std::string("unknown residue: ") + res_name;
		}

		inline int res_mol_type(const S &res_name) {
			int i = 0;
			for (auto && mol_type : { "RNA", "DNA", "protein" }) {
				const Names &names = Names::instance(mol_type);
				auto it = std::find_if(names.res.begin(), names.res.end(), [&res_name, &names](const S &s) {
					const names_t &v = names.alias.at(s);
					return s == res_name || std::find(v.begin(), v.end(), res_name) != v.end();
				});
				if (it != names.res.end()) {
					return i;
				}
				i++;
			}
			throw std::string("unknown residue: ") + res_name;
		}

		template<typename T>
		inline T res_mol_type(const S &res_name) {
			for (auto && mol_type : { "RNA", "DNA", "protein" }) {
				//std::cout << mol_type << std::endl;
				const Names &names = Names::instance(mol_type);
				auto it = std::find_if(names.res.begin(), names.res.end(), [&res_name, &names](const S &s) {
					const names_t &v = names.alias.at(s);
					return s == res_name || std::find(v.begin(), v.end(), res_name) != v.end();
				});
				if (it != names.res.end()) {
					return mol_type;
				}
			}
			throw std::string("unknown residue: ") + res_name;
		}

		inline const names_t &res_included_atoms(const S &res_name) {
			return Names::instance(res_mol_type<std::string>(res_name)).atoms_res.at(res_name);
		}

	}
END_JN

