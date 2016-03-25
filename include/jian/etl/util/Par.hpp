#ifndef PAR_H
#define PAR_H

#include "../std.hpp"

namespace jian {

class Par {
public:
    std::map<std::string, std::deque<std::string>> _pars;
    std::deque<std::string> _orig_pars;

    Par() {}

    Par(int argc, char **argv) {
        read(argc, argv);
    }

    Par(std::string str) {
        read(str);
    }

    void read(int argc, char **argv) {
        for (int i = 0; i < argc; i++) _orig_pars.push_back(argv[i]);
        std::string key;
        std::deque<std::string> values;
        int n = 0;
        for (int i = 1; i < argc; i++) {
            if (argv[i][0] == '-') {
                if (n != 0) {
                    _pars[key] = values;
                } else {
                    _pars["global"] = values;
                }
                std::string str(argv[i]);
                key = str.substr(1, str.size() - 1);
                values.clear();
                n++;
            } else {
                values.push_back(argv[i]);
            }
        }
        if (n != 0) {
            _pars[key] = values;
        } else {
            _pars["global"] = values;
        }

    }

    std::deque<std::string> &operator [](const std::string &s) {
        return _pars.at(s);
    }

    const std::deque<std::string> &operator [](const std::string &s) const {
        return _pars.at(s);
    }

    std::string &operator [](int n) {
        return _orig_pars[n];
    }

    const std::string &operator [](int n) const {
        return _orig_pars[n];
    }

    bool count(const std::string &s) const {
        return _pars.find(s) != _pars.end();
    }

    bool has(const std::string &s) const {
        return _pars.find(s) != _pars.end();
    }

    void read(string par_file) {
        std::ifstream ifile(par_file.c_str());
        if (!ifile) throw "Par::read(string) error! Can't open file '" + par_file + "'.";
        std::string key, value;
        while (ifile >> key >> value) {
            _pars[key].push_back(value);
        }
        ifile.close();
    }

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
        return boost::lexical_cast<T>(s);
    }

    template<typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 42>
    T parse(const std::string &s) const {
        return boost::lexical_cast<T>(s);
    }

    template<typename T, std::enable_if_t<std::is_same<std::decay_t<T>, std::string>::value, int> = 42>
    T parse(const std::string &s) const {
        return s;
    }

};

} /// namespace jian

#endif
