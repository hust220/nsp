#pragma once

#include <string>
#include <fstream>
#include <deque>
#include <map>
#include <list>
#include <vector>
#include <iostream>
#include "string.hpp"
#include "traits.hpp"

namespace jian {

class Par {
public:
    using pars_t = std::deque<std::string>;

    std::map<std::string, pars_t> _pars;
    pars_t _orig_pars;

	Par() = default;

	Par(const Par &par) = default;

    Par(int argc, char **argv);

    Par(str_t str);

	template<typename _FirstVal, typename... _RestVals>
	Par(const std::string &key, _FirstVal &&first_val, _RestVals &&...rest_vals) {
		(*this)(key, first_val, rest_vals...);
	}

    void read(int argc, char **argv);

    std::deque<std::string> &operator [](const std::string &s);

    const std::deque<std::string> &operator [](const std::string &s) const;

    std::deque<std::string> &operator [](const std::vector<std::string> &s);

    const std::deque<std::string> &operator [](const std::vector<std::string> &s) const;

    std::string &operator [](int n);

    const std::string &operator [](int n) const;

	Par &operator ()(const std::string &key) {
		_pars[key];
		return *this;
	}

	template<typename _FirstVal, typename... _RestVals>
	Par &operator ()(const std::string &key, _FirstVal &&first_val, _RestVals &&...rest_vals) {
		_pars[key].push_back(lexical_cast<std::string>(first_val));
		return (*this)(key, rest_vals...);
	}

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

	std::list<str_t> &keys_chain() const {
		static std::list<str_t> chain;
		return chain;
	}

    std::string get() const {
		std::ostringstream stream;
		stream << "jian::Par::get error! Didn't found parameters for keys:";
		for (auto && key : keys_chain()) stream << " " << key;
		throw stream.str();
	}

    template<typename K, typename... V>
    std::string get(K &&s, V && ...pars) const {
        if (_pars.count(s)) {
			keys_chain().clear();
            return _pars.at(s)[0];
        } else {
			keys_chain().push_back(s);
            return get(pars...);
        }
    }

	pars_t getv() const {
		std::ostringstream stream;
		stream << "jian::Par::getv error! Didn't found parameters for keys:";
		for (auto && key : keys_chain()) stream << " " << key;
		throw stream.str();
	}

	template<typename K, typename... V>
	pars_t getv(K &&s, V && ...pars) const {
		if (_pars.count(s)) {
			keys_chain().clear();
			return _pars.at(s);
		}
		else {
			keys_chain().push_back(s);
			return getv(pars...);
		}
	}

	template<typename T>
	T parse(const std::string &s) const {
		return jian::lexical_cast<T>(s);
	}

};

} // namespace jian

