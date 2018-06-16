#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <algorithm>
#include "file.hpp"
#include "log.hpp"
#include "pdb_parser.hpp"
#include "pdb_residue.hpp"
#include "pdb_names.hpp"
#include "pdb_io.hpp"

namespace jian {

bool is_mol_type(const Residue &res, Str type = "") {
	if (type == "") {
		return true;
	}
	else {
		const pdb::Names & v = pdb::Names::instance(type);
		S name = jian::upper(res.name);

		return std::find_if(v.res.begin(), v.res.end(), [&name, &v](const S & s) {
			const pdb::names_t & l = v.alias.at(s);
			return name == s || std::find(l.begin(), l.end(), name) != l.end();
		}) != v.res.end();
	}
}

auto get_sort_keys() {
	Map<std::string, std::vector<std::string>> keys{
		{"A",{"P","O1P","O2P",
			  "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
			  "N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"}},
		{"U",{"P","O1P","O2P",
			  "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
			  "N1","C2","O2","N3","C4","O4","C5","C6"}},
		{"G",{"P","O1P","O2P",
			  "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
			  "N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"}},
		{"C",{"P","O1P","O2P",
			  "O5*","C5*","C4*","O4*","C3*","O3*","C2*","O2*","C1*",
			  "N1","C2","O2","N3","C4","N4","C5","C6"}}
	};
	std::map<std::string, std::map<std::string, int>> sort_keys;
	int index = 0;
	for (auto && res_name : { "A", "U", "G", "C" }) {
		for (int i = 0; i < keys[res_name].size(); i++) {
			sort_keys[res_name][keys[res_name][i]] = index;
			index++;
		}
	}
	return sort_keys;
}

void sort(Residue &res) {
	static auto sort_keys = get_sort_keys();
	auto & keys = sort_keys;
	std::sort(res.begin(), res.end(), [&res, &keys](auto &&a1, auto &&a2) {
		return keys[res.name][a1.name] < keys[res.name][a2.name];
	});
}

Str format_res_name(const Str &s) {
	std::smatch result;
	if (std::regex_match(s, result, std::regex("^(\\w+)\\d+$"))) {
		// tLeap would append a '5' after the name of first residue
		// and '3' after the name of the last residue
		return result[1];
	}
	else return s;
}

Atom &Residue::operator [](int n) {
	return std::deque<Atom>::operator [](n);
}

const Atom &Residue::operator [](int n) const {
	return std::deque<Atom>::operator [](n);
}

Atom &Residue::operator [](const S &s) {
	for (auto &&atom : *this) if (atom.name == s) { return atom; }
	std::cerr << "jian::Residue::operator[] error! Not found atom '" << s << "'!" << std::endl;
	std::cerr << *this << std::endl;
	die();
}

const Atom &Residue::operator[](const S &s) const {
	for (auto &&atom : *this) if (atom.name == s) { return atom; }
	std::cerr << "jian::Residue::operator[] error! Not found atom '" << s << "'!" << std::endl;
	std::cerr << *this << std::endl;
	die();
}

Residue &Residue::operator+=(const Residue &res) {
	for (auto && atom : res) {
		push_back(atom);
	}
	return *this;
}

}

