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

BEGIN_JN

class Par {
public:
    using pars_t = std::deque<Str>;

    std::map<Str, pars_t> _pars;
    pars_t _orig_pars;

    Int m_argc;
    Char **m_argv;

	Par() = default;

	Par(const Par &par) = default;

    Par(int argc, char **argv);

    Par(Str str);

	template<typename _FirstVal, typename... _RestVals>
	Par(const Str &key, _FirstVal &&first_val, _RestVals &&...rest_vals) {
		(*this)(key, first_val, rest_vals...);
	}

    void read(int argc, char **argv);

    std::deque<Str> &operator [](const Str &s);

    const std::deque<Str> &operator [](const Str &s) const;

    std::deque<Str> &operator [](const std::vector<Str> &s);

    const std::deque<Str> &operator [](const std::vector<Str> &s) const;

    Str &operator [](int n);

    const Str &operator [](int n) const;

	Par &operator ()(const Str &key) {
		_pars[key];
		return *this;
	}

	template<typename _FirstVal, typename... _RestVals>
	Par &operator ()(const Str &key, _FirstVal &&first_val, _RestVals &&...rest_vals) {
		_pars[key].push_back(lexical_cast<Str>(first_val));
		return (*this)(key, rest_vals...);
	}

	bool count(const Str &s) const;

    void read(Str par_file);

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

	template<typename T, typename _First, typename... _Rest>
	void setv(T &&t, _First &&first, _Rest && ...rest) const {
		if (_pars.find(first) != _pars.end()) {
			auto && l = _pars.at(first);
			t.resize(l.size());
			std::copy(l.begin(), l.end(), t.begin());
		}
		else {
			setv(t, rest...);
		}
	}

	std::list<Str> &keys_chain() const {
		static std::list<Str> chain;
		return chain;
	}

    Str get() const {
		std::ostringstream stream;
		stream << "jian::Par::get error! Didn't found parameters for keys:";
		for (auto && key : keys_chain()) stream << " " << key;
		throw stream.str();
	}

    template<typename K, typename... V>
    Str get(K &&s, V && ...pars) const {
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
	T parse(const Str &s) const {
		return jian::lexical_cast<T>(s);
	}

};

END_JN

