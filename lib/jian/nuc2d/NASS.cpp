#include <numeric>
#include <deque>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <algorithm>
#include <mutex>
#include "../pp.hpp"
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "NASS.hpp"

BEGIN_JN

	namespace nass_detail {

		std::map<char, char> get_map_keys() {
			std::map<char, char> map_keys;
			for (auto &&pair : NASS::instance().paired_keys) {
				map_keys[pair.first] = pair.second;
				map_keys[pair.second] = pair.first;
			}
			return map_keys;
		}

		std::map<char, int> get_pos_keys() {
			std::map<char, int> pos_keys;
			int index = 1;
			for (auto &&pair : NASS::instance().paired_keys) {
				pos_keys[pair.first] = index;
				pos_keys[pair.second] = -index;
				index++;
			}
			return pos_keys;
		}

		std::map<char, char> map_keys = get_map_keys();

		std::map<char, int> pos_keys = get_pos_keys();

		std::mutex mt;

	} // namespace nass_detail

	NASS::NASS() {
		Str name = Env::lib() + "/RNA/pars/nuc2d/dbn.symbols";
		FileLines lines(name);
		for (auto it = lines.begin(); it != lines.end(); it++) {
			if (it->arr.size() > 1) {
				if (it->arr[1] == "paired_keys") {
					it++;
					for (int i = 0; i < it->arr.size() / 2; i++) {
						paired_keys.push_back({ it->arr[2 * i][0], it->arr[2 * i + 1][0] });
					}
				}
				else if (it->arr[1] == "unpaired_keys") {
					it++;
					for (auto && s : it->arr) {
						unpaired_keys.push_back(s[0]);
					}
				}
				else if (it->arr[1] == "break_keys") {
					it++;
					for (auto && s : it->arr) {
						break_keys.push_back(s[0]);
					}
				}
			}
		}
	}

	std::vector<int> NASS::get_bps(const Str &ss) {
		const NASS &nass = NASS::instance();
		std::deque<std::deque<int>> dq;
		std::vector<int> bps;
		int i, j, k;

		dq.resize(nass.paired_keys.size());
		bps.resize(std::accumulate(ss.begin(), ss.end(), 0, [&nass](size_t n, char c)->size_t {
			return n + (std::find(nass.break_keys.begin(), nass.break_keys.end(), c)!= nass.break_keys.end() ? 0 : 1);
		}), -1);
		i = 0;
		for (auto && c : ss) {
			if (std::find(nass.break_keys.begin(), nass.break_keys.end(), c)!= nass.break_keys.end()) {
				continue;
			}
			auto it1 = std::find_if(nass.paired_keys.begin(), nass.paired_keys.end(), [c](auto && p){
				return p.first == c;
			});
			auto it2 = std::find_if(nass.paired_keys.begin(), nass.paired_keys.end(), [c](auto && p) {
				return p.second == c;
			});
			if (it1 != nass.paired_keys.end()) {
				j = std::distance(nass.paired_keys.begin(), it1);
				dq[j].push_back(i);
			}
			else if (it2 != nass.paired_keys.end()) {
				j = std::distance(nass.paired_keys.begin(), it2);
				k = dq[j].back();
				bps[k] = i;
				bps[i] = k;
				dq[j].pop_back();
			}
			i++;
		}
		return bps;
	}

	const NASS &NASS::instance() {
		static NASS nuc_ss;
		return nuc_ss;
	}

	bool NASS::is_char_ss(char c) {
		return
			std::find_if( instance().paired_keys.begin(), instance().paired_keys.end(), [&](const std::pair<char, char> &pair) {
				return pair.first == c || pair.second == c;
			}) != instance().paired_keys.end() ||
			std::find(instance().unpaired_keys.begin(), instance().unpaired_keys.end(), c) != instance().unpaired_keys.end() ||
			std::find(instance().break_keys.begin(), instance().break_keys.end(), c) != instance().break_keys.end();
	}


	bool NASS::check_ss(const Str &ss) {
		Str info_errors;
		return check_ss(ss, info_errors);
	}

	bool NASS::check_ss(const Str &ss, Str &info_errors) {
		std::lock_guard<std::mutex> gd(nass_detail::mt);
		std::ostringstream stream;

		for (auto &&s : ss) {
			if (!is_char_ss(s)) {
				stream << "Illegal character in secondary structure: '" << s << "'.\n";
				info_errors = stream.str();
				return false;
			}
		}

		std::map<char, int> map;
		for (auto &&s : ss) {
			if (nass_detail::pos_keys.count(s)) {
				if (nass_detail::pos_keys[s] > 0) {
					if (map.count(s)) map[s]++; else map[s] = 1;
				}
				else {
					char c = nass_detail::map_keys[s];
					if (map.count(c)) {
						map[c]--;
						if (map[c] < 0) return false;
					}
					else {
						return false;
					}
				}
			}
		}

		for (auto &&pair : map) {
			if (pair.second != 0) {
				info_errors = 
					"The number of left parenthesis (or squares) should equal to the number of right parenthesis (or squares)!\n";
				return false;
			}
		}

		return true;
	}

	int NASS::len_ss(const Str &ss) {
		std::lock_guard<std::mutex> gd(nass_detail::mt);

		return std::count_if(ss.begin(), ss.end(), [&](const char &c) {
			return std::count_if(instance().paired_keys.begin(), instance().paired_keys.end(), [&](const std::pair<char, char> &pair) {
				return pair.first == c || pair.second == c;
			}) || std::count(instance().unpaired_keys.begin(), instance().unpaired_keys.end(), c);
		});
	}

	Str NASS::pure_ss(const Str &ss) {
		std::lock_guard<std::mutex> gd(nass_detail::mt);

		Str p_ss;
		std::copy_if(ss.begin(), ss.end(), std::back_inserter(p_ss), [](const char &c) {return c != '&'; });
		return p_ss;
	}

	Str NASS::lower_ss(const Str &ss, int n) {
		std::lock_guard<std::mutex> gd(nass_detail::mt);

		Str l_ss;
		std::transform(ss.begin(), ss.end(), std::back_inserter(l_ss), [&](const char &c) {
			if (nass_detail::pos_keys.count(c) && (nass_detail::pos_keys[c] > n || nass_detail::pos_keys[c] < -n)) {
				return '.';
			}
			else {
				return c;
			}
		});
		return l_ss;
	}

	Str NASS::hinge_ss(const Str &ss) {
		Str h_ss;
		for (auto && i : ss) {
			if (i == '(' || i == ')') {
				h_ss += i;
			}
		}
		//EACH(i, ss, if (i == '(' || i == ')') h_ss += i);
		return h_ss;
	}

	bool NASS::seq_match_ss(const Str &seq, const Str &ss) {
		return std::count_if(seq.begin(), seq.end(), [](auto && c) {return c != '&'; }) ==
			std::count_if(ss.begin(), ss.end(), [](auto && c) {return c != '&'; });
	}

END_JN

