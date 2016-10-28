#pragma once

#include <string>
#include <fstream>
#include <deque>
#include <map>
#include <vector>
#include <iostream>
#include "string.hpp"

namespace jian {

class Par {
public:
    using pars_t = std::deque<std::string>;

    std::map<std::string, pars_t> _pars;
    pars_t _orig_pars;

    Par() {}
    Par(int argc, char **argv);
    Par(std::string str);
    void read(int argc, char **argv);
    std::deque<std::string> &operator [](const std::string &s);
    const std::deque<std::string> &operator [](const std::string &s) const;
    std::deque<std::string> &operator [](const std::vector<std::string> &s);
    const std::deque<std::string> &operator [](const std::vector<std::string> &s) const;
    std::string &operator [](int n);
    const std::string &operator [](int n) const;
    bool count(const std::string &s) const;
    void read(std::string par_file);
    friend std::ostream &operator <<(std::ostream &out, const Par &par);

	template<typename K>
	bool has(K &&k) const {
		return _pars.find(k) != _pars.end();
	}

	template<typename K, typename U, typename... V>
	bool has(K &&k, U &&u, V && ...rest) const {
		return has(k) || has(u, rest...);
	}

	template<typename T>
	void set(T &&v) const {}

	template<typename T, typename K, typename... V>
	void set(T &&v, K &&s, V && ...rest) const {
		if (_pars.find(s) != _pars.end()) {
			v = parse<std::decay_t<T>>(_pars.at(s)[0]);
		}
		else set(v, rest...);
	}

	template<typename T>
	void setv(T &&v) const {}

	template<typename T, typename K, typename... V>
	void setv(T &&t, K &&k, V && ...rest) const {
		if (_pars.find(k) != _pars.end()) {
			auto && l = _pars.at(k);
			t.resize(l.size());
			std::copy(l.begin(), l.end(), t.begin());
		}
		else {
			setv(t, rest...);
		}
	}

    std::string get() const {
        throw "jian::Par::get error! Didn't found parameters for keys!";
    }

    template<typename K, typename... V>
    std::string get(K &&s, V && ...pars) const {
        if (_pars.count(s)) {
            return _pars.at(s)[0];
        } else {
            return get(pars...);
        }
    }

	pars_t getv() const {
		throw "jian::Par::getv error! Didn't found parameters for keys!";
	}

	template<typename K, typename... V>
	pars_t getv(K &&s, V && ...pars) const {
		if (_pars.count(s)) {
			return _pars.at(s);
		}
		else {
			return getv(pars...);
		}
	}

	template<typename T>
	T parse(const std::string &s) const {
		return jian::lexical_cast<T>(s);
	}

};

} // namespace jian

