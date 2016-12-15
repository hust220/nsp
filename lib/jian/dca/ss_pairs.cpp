#include <array>
#include <list>
#include <string>
#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include "../utils/file.hpp"
#include "../nuc2d.hpp"
#include "ss_pairs.hpp"
#include "../utils/file.hpp"

BEGIN_JN
namespace dca {

	std::vector<std::array<char, 2>> patterns = { {'.', '.'}, {'(', ')'}, {'[', ']'}, {'{', '}'} };

void pairs_sort(pairs_t &pairs) {
    pairs.sort([](auto && pair1, auto && pair2){
		return pair1[0] < pair2[0] || (pair1[0] == pair2[0] && pair1[1] < pair2[1]);
	});
}

pairs_t pairs_from_file(const Str &file_name, int size) {
    auto && tuples = tuples_from_file(file_name, size);
    pairs_t pairs;
    for (auto && tuple : tuples) {
        pairs.push_back({tuple.a, tuple.b});
    }
    return pairs;
}

tuples_t tuples_from_file(const Str &file_name, int size) {
    tuples_t tuples;
	int a, b;
	Num c;

	for (auto &&l : FileLines(file_name)) {
		if (JN_ size(l.arr) == 2) {
			a = JN_INT(l.arr[0]) - 1;
			b = JN_INT(l.arr[1]) - 1;
			c = 1;
		}
		else if (JN_ size(l.arr) >= 3) {
			a = JN_INT(l.arr[0]) - 1;
			b = JN_INT(l.arr[1]) - 1;
			c = JN_NUM(l.arr[2]);
		}
		if (a >= 0 && b - a - 1 >= 4) {
			auto p = std::minmax(a, b);
			if (c < tuples.back().c) {
				if (tuples.size() >= size) {
				}
				else {
					tuples.push_back({ p.first, p.second, c });
				}
			}
			else if (c > tuples.front().c) {
				tuples.push_front({ p.first, p.second, c });
			}
			else {
				auto it_p = tuples.begin();
				for (auto it = std::next(it_p); it != tuples.end(); it = std::next(it)) {
					if (c > it->c && c < it_p->c) {
						tuples.insert(it, { p.first, p.second, c });
						break;
					}
					it_p = it;
				}
			}
			if (tuples.size() > size) {
				tuples.pop_back();
			}
		}
	}
    return tuples;
}

pairs_t pairs_from_ss(const ss_t &ss) {
    pairs_t pairs;
    auto & keys = NASS::instance().paired_keys;
    std::vector<std::deque<int>> v(keys.size());
    int i = 0;
    for (auto && c : ss) {
        if (c == '&') continue;
        auto it1 = std::find_if(keys.begin(), keys.end(), [&c](auto && key){return key.first == c;});
        if (it1 != keys.end()) {
            int n = std::distance(keys.begin(), it1);
            v[n].push_back(i);
        } else {
            auto it2 = std::find_if(keys.begin(), keys.end(), [&c](auto && key){return key.second == c;});
            if (it2 != keys.end()) {
                int n = std::distance(keys.begin(), it2);
                pairs.push_back({v[n].back(), i});
                v[n].pop_back();
            } else {
            }
        }
        i++;
    }
    return pairs;
}

void print_pairs(const pairs_t & pairs) {
    for (auto && pair : pairs) {
        std::cout << pair[0]+1 << ' ' << pair[1]+1 << std::endl;
    }
}

void print_tuples(const tuples_t & tuples) {
    for (auto && tuple : tuples) {
        std::cout << tuple.a+1 << ' ' << tuple.b+1 << ' ' << tuple.c << std::endl;
    }
}

int get_level(const ss_t &ss, const pair_t &l) {
    if (ss[l[0]] != '.' || ss[l[1]] != '.') {
        return 0;
    }
    for (int level = 1; level < patterns.size(); level++) {
        int score = 0, flag = 1;
        for (int i = l[0]+1; i < l[1]; i++) {
            if (ss[i] == patterns[level][0]) {
                score++;
            } else if (ss[i] == patterns[level][1]) {
                score--;
            }
            if (score < 0) {
                flag = 0;
                break;
            }
        }
        if (score != 0) {
            flag = 0;
        }
        if (flag == 1) {
            return level;
        }
    }
    return 0;
}

ss_t pairs_to_ss(const pairs_t &pairs, int size) {
    ss_t ss(size, '.');
    int level;
    for (auto && l : pairs) {
        level = get_level(ss, l);
        if (level != 0) {
            ss[l[0]] = patterns[level][0];
            ss[l[1]] = patterns[level][1];
        }
    }
    return ss;
}


} // namespace dca
END_JN

