#pragma once

#include <string>
#include <fstream>
#include <deque>
#include <map>
#include <vector>

namespace jian {

class Par {
public:
    using val_t = std::deque<std::string>;

    std::map<std::string, val_t> _pars;
    val_t _orig_pars;

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
    bool has(const std::string &s) const;
    void read(std::string par_file);
    friend std::ostream &operator <<(std::ostream &out, const Par &par);

    template<typename T>
    void set(T &&v) const {}

    template<typename T, typename K, typename... V>
    void set(T &&v, K &&s, V && ...pars) const {
        if (_pars.count(s)) {
            v = parse<std::decay_t<T>>(_pars.at(s)[0]);
        } else set(v, pars...);
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

