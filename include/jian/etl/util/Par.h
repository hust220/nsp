#ifndef PAR_H
#define PAR_H

#include "../std.h"

namespace jian {

class Par: public std::map<std::string, std::vector<std::string>> {
public:
    Par() {}

    Par(int argc, char **argv) {
        read(argc, argv);
    }

    Par(string str) {
        read(str);
    }

    void read(int argc, char **argv) {
        std::string key;
        std::vector<std::string> values;
        map<string, vector<string>> par;
        int n = 0;
        for (int i = 1; i < argc; i++) {
            if (argv[i][0] == '-') {
                if (n != 0) {
                    (*this)[key] = values;
                } else {
                    (*this)["global"] = values;
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
            (*this)[key] = values;
        } else {
            (*this)["global"] = values;
        }

    }

    void read(string par_file) {
        map<string, vector<string>> par;
        ifstream ifile(par_file.c_str());
        if (!ifile) throw "Par::read(string) error! Can't open file '" + par_file + "'.";
        string key, value;
        while (ifile >> key >> value) {
            (*this)[key].push_back(value);
        }
        ifile.close();
    }

};

} /// namespace jian

#endif
