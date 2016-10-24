#pragma once

#include <string>
#include <fstream>
#include <deque>
#include <map>
#include <vector>
#include <iostream>

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

    bool has(const std::string &s) const;

    template<typename T>
    void set(T &&v) const {}

    template<typename T, typename K, typename... V>
    void set(T &&v, K &&s, V && ...pars) const {
        if (_pars.count(s)) {
            v = parse<std::decay_t<T>>(_pars.at(s)[0]);
        } else set(v, pars...);
    }

    pars_t getall() const {
        throw "jian::Par::get error! Didn't found parameters for keys!";
    }

    template<typename K, typename... V>
    pars_t getall(K &&s, V && ...pars) const {
        if (_pars.count(s)) {
            return _pars.at(s);
        } else {
            return getall(pars...);
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

    template<typename T, std::enable_if_t<std::is_integral<T>::value, int> = 42>
    T parse(const std::string &s) const {
        return std::stoi(s);
    }

    template<typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 42>
    T parse(const std::string &s) const {
        return std::stod(s);
    }

    template<typename T, std::enable_if_t<std::is_same<std::decay_t<T>, std::string>::value, int> = 42>
    T parse(const std::string &s) const {
        return s;
    }

};

} // namespace jian

