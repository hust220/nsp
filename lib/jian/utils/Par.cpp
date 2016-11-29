#include <iostream>
#include <cctype>

#include "Par.hpp"
#include "file.hpp"

namespace jian {

Par::Par(int argc, char **argv) {
    read(argc, argv);
}

Par::Par(str_t str) {
    read(str);
}

void Par::read(str_t par_file) {
	BEGIN_READ_FILE(par_file, " ") {
		if (!(F.empty())) {
			_pars[F[0]] = pars_t();
			for (auto it = std::next(F.begin()); it != F.end(); it++) {
				_pars[F[0]].push_back(*it);
			}
			if (F[0] == "par") {
				for (auto && p : _pars["par"]) {
					read(p);
				}
			}
		}
	} END_READ_FILE;
}

void Par::read(int argc, char **argv) {
    for (int i = 0; i < argc; i++) _orig_pars.push_back(argv[i]);
    str_t key;
    std::deque<str_t> values;
    int n = 0;
    for (int i = 1; i < argc; i++) {
        if (str_t(argv[i]).size() >=2 && argv[i][0] == '-' && std::isalpha(argv[i][1])) {
            if (n != 0) {
                _pars[key] = values;
                if (key == "par") {
                    for (auto && p : values) {
                        read(p);
                    }
                }
            } else {
                _pars["global"] = values;
            }
            str_t str(argv[i]);
            key = str.substr(1, str.size() - 1);
            values.clear();
            n++;
        } else {
            values.push_back(argv[i]);
        }
    }
    if (n != 0) {
        _pars[key] = values;
        if (key == "par") {
            for (auto && p : values) {
                read(p);
            }
        }
    } else {
        _pars["global"] = values;
    }
//    if (has("par")) {
//        for (auto && p : (*this)["par"]) {
//            read(p);
//        }
//    }

}

std::deque<str_t> &Par::operator [](const str_t &s) {
    return _pars.at(s);
}

const std::deque<str_t> &Par::operator [](const str_t &s) const {
    return _pars.at(s);
}

std::deque<str_t> &Par::operator [](const std::vector<str_t> &v) {
    for (auto && s : v) {
        if (has(s)) {
            return _pars.at(s);
        }
    }
    throw "jian::Par error!";
}

const std::deque<str_t> &Par::operator [](const std::vector<str_t> &v) const {
    for (auto && s : v) {
        if (has(s)) {
            return _pars.at(s);
        }
    }
    throw "jian::Par error!";
}

str_t &Par::operator [](int n) {
    return _orig_pars[n];
}

const str_t &Par::operator [](int n) const {
    return _orig_pars[n];
}

bool Par::count(const str_t &s) const {
    return _pars.find(s) != _pars.end();
}

std::ostream &operator <<(std::ostream &out, const Par &par) {
    for (auto && i : par._pars) {
        out << i.first;
        for (auto && j : i.second) {
            out << ' ' << j;
        }
        out << std::endl;
    }
    return out;
}

} // namespace jian

