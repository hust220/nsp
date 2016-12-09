#include <iostream>
#include <cctype>

#include "Par.hpp"
#include "file.hpp"

BEGIN_JN

Par::Par(int argc, char **argv) {
    read(argc, argv);
}

Par::Par(Str str) {
    read(str);
}

void Par::read(Str par_file) {
	for (auto &&it : FileLines(par_file)) {
		if (!(it.arr.empty())) {
			_pars[it.arr[0]] = pars_t();
			for (auto it2 = std::next(it.arr.begin()); it2 != it.arr.end(); it2++) {
				_pars[it.arr[0]].push_back(*it2);
			}
			if (it.arr[0] == "par") {
				for (auto && p : _pars["par"]) {
					read(p);
				}
			}
		}
	}
}

void Par::read(int argc, char **argv) {
    for (int i = 0; i < argc; i++) _orig_pars.push_back(argv[i]);
    Str key;
    std::deque<Str> values;
    int n = 0;
    for (int i = 1; i < argc; i++) {
        if (Str(argv[i]).size() >=2 && argv[i][0] == '-' && std::isalpha(argv[i][1])) {
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
            Str str(argv[i]);
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

std::deque<Str> &Par::operator [](const Str &s) {
    return _pars.at(s);
}

const std::deque<Str> &Par::operator [](const Str &s) const {
    return _pars.at(s);
}

std::deque<Str> &Par::operator [](const std::vector<Str> &v) {
    for (auto && s : v) {
        if (has(s)) {
            return _pars.at(s);
        }
    }
    throw "jian::Par error!";
}

const std::deque<Str> &Par::operator [](const std::vector<Str> &v) const {
    for (auto && s : v) {
        if (has(s)) {
            return _pars.at(s);
        }
    }
    throw "jian::Par error!";
}

Str &Par::operator [](int n) {
    return _orig_pars[n];
}

const Str &Par::operator [](int n) const {
    return _orig_pars[n];
}

bool Par::count(const Str &s) const {
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

END_JN

