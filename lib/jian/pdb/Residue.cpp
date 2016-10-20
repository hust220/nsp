#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <algorithm>
#include "../utils/file.hpp"
#include "molstream.hpp"
#include "Residue.hpp"

namespace jian {

	bool res_is_type(const Residue &res, std::string type) {
		static std::vector<std::string> v{ "A", "U", "G", "C",
										   "RA", "RU", "RG", "RC",
										   "A5", "U5", "G5", "C5",
										   "A3", "U3", "G3", "C3" };
		if (type == "RNA") {
			if (std::find(v.begin(), v.end(), jian::upper(res.name)) != v.end()) {
				return true;
			}
			else {
				return false;
			}
		}
		else {
			return true;
		}
	}

	auto Residue::get_sort_keys() {
		thread_local static std::map<std::string, std::vector<std::string>> keys{
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

	void Residue::sort() {
		thread_local static auto sort_keys = get_sort_keys();
		auto & keys = sort_keys;
		std::sort(this->begin(), this->end(), [&](auto &&a1, auto &&a2) {
			return keys[name][a1.name] < keys[name][a2.name];
		});
	}


	std::string Residue::format_name(const std::string &s) {
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

	Atom &Residue::operator [](const std::string &s) {
		for (auto &&atom : *this) if (atom.name == s) { return atom; }
		throw "jian::Residue::operator[] error! Not found atom!";
	}

	const Atom &Residue::operator [](const std::string &s) const {
		for (auto &&atom : *this) if (atom.name == s) { return atom; }
		throw "jian::Residue::operator[] error! Not found atom!";
	}

} // namespace jian

