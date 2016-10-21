#include <vector>
#include <utility>
#include <map>
#include <string>
#include <algorithm>
#include <mutex>
#include "../pp.hpp"
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "NucSS.hpp"

namespace jian {

namespace nass_detail {

std::map<char, char> get_map_keys() {
    std::map<char, char> map_keys;
    for (auto &&pair: NucSS::instance().paired_keys) {
        map_keys[pair.first] = pair.second;
        map_keys[pair.second] = pair.first;
    }
    return map_keys;
}

std::map<char, int> get_pos_keys() {
    std::map<char, int> pos_keys;
    int index = 1;
    for (auto &&pair: NucSS::instance().paired_keys) {
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

NucSS::NucSS() {
    std::string name = Env::lib() + "/RNA/pars/nuc2d/dbn.symbols";
	EACH_SPLIT_LINE(name.c_str(), " ",
		if (F.size() > 1) {
			if (F[1] == "paired_keys") {
				std::getline(ifile, L);
				jian::tokenize(L, F, " ");
				for (int i = 0; i < F.size() / 2; i++) {
					paired_keys.push_back({ F[2 * i][0], F[2 * i + 1][0] });
				}
			} else if (F[1] == "unpaired_keys") {
				std::getline(ifile, L);
				jian::tokenize(L, F, " ");
				for (auto && s : F) {
					unpaired_keys.push_back(s[0]);
				}
			} else if (F[1] == "break_keys") {
				std::getline(ifile, L);
				jian::tokenize(L, F, " ");
				for (auto && s : F) {
					break_keys.push_back(s[0]);
				}
			}
		}
    );
}

NucSS &NucSS::instance() {
    static NucSS nuc_ss;
    return nuc_ss;
}

bool NucSS::check_ss(const std::string &ss) {
    std::lock_guard<std::mutex> gd(nass_detail::mt);

    for (auto &&s: ss) {
        if (! std::count_if(instance().paired_keys.begin(), instance().paired_keys.end(), [&](const std::pair<char, char> &pair) {
                return pair.first == s || pair.second == s;
            }) && ! std::count(instance().unpaired_keys.begin(), instance().unpaired_keys.end(), s) 
               && ! std::count(instance().break_keys.begin(), instance().break_keys.end(), s)) {
            return false;
        }
    }
    std::map<char, int> map;
    for (auto &&s: ss) {
        if (nass_detail::pos_keys.count(s)) {
            if (nass_detail::pos_keys[s] > 0) {
                if (map.count(s)) map[s]++; else map[s] = 1;
            } else {
                char c = nass_detail::map_keys[s];
                if (map.count(c)) {
                    map[c]--;    
                    if (map[c] < 0) return false;
                } else {
                    return false;    
                }
            }
        }
    }

    for (auto &&pair: map) if (pair.second != 0) return false;

    return true;
}

int NucSS::len_ss(const std::string &ss) {
    std::lock_guard<std::mutex> gd(nass_detail::mt);

    return std::count_if(ss.begin(), ss.end(), [&](const char &c){
        return std::count_if(instance().paired_keys.begin(), instance().paired_keys.end(), [&](const std::pair<char, char> &pair){
            return pair.first == c || pair.second == c;
        }) || std::count(instance().unpaired_keys.begin(), instance().unpaired_keys.end(), c);
    });
}

std::string NucSS::pure_ss(const std::string &ss) {
    std::lock_guard<std::mutex> gd(nass_detail::mt);

    std::string p_ss;
    std::copy_if(ss.begin(), ss.end(), std::back_inserter(p_ss), [](const char &c){return c != '&';});
    return p_ss;
}
    
std::string NucSS::lower_ss(const std::string &ss, int n) {
    std::lock_guard<std::mutex> gd(nass_detail::mt);

    std::string l_ss;
    std::transform(ss.begin(), ss.end(), std::back_inserter(l_ss), [&](const char &c){
        if (nass_detail::pos_keys.count(c) && (nass_detail::pos_keys[c] > n || nass_detail::pos_keys[c] < -n)) {
            return '.';    
        } else {
            return c;
        }
    });
    return l_ss;
}

std::string NucSS::hinge_ss(const std::string &ss) {
    std::string h_ss; EACH(i, ss, if (i == '(' || i == ')') h_ss += i); return h_ss;
}

bool NucSS::seq_match_ss(const std::string &seq, const std::string &ss) {
	return std::count_if(seq.begin(), seq.end(), [](auto && c) {return c != '&'; }) ==
		std::count_if(ss.begin(), ss.end(), [](auto && c) {return c != '&'; });
}

} // namespace jian

