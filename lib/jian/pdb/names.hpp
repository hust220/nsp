#pragma once

#include <string>
#include <vector>
#include <map>

namespace jian {
	namespace pdb {

		using names_t = std::vector<std::string>;
		using map_names_t = std::map<std::string, names_t>;

		int res_type(const std::string &res_name);

		class Names {
		public:
			static const Names & instance(std::string mol_type);

			names_t res;
			map_names_t alias;
			names_t atoms_phos;
			names_t atoms_sugar;
			map_names_t atoms_base;
			map_names_t atoms_res;
		};

	}
} // namespace jian

