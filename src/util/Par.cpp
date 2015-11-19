#include "Par.h"

namespace jian {

//map<string, vector<string>> Par::read(int argc, char **argv) {
void Par::read(int argc, char **argv) {
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

//    return par;
}

void Par::read(string par_file) {
//map<string, vector<string>> Par::read(string par_file) {
    map<string, vector<string>> par;
    ifstream ifile(par_file.c_str());
    ifile || die("Par::read(string) error! Can't open file '" + par_file + "'.");
    string key, value;
    while (ifile >> key >> value) {
        (*this)[key].push_back(value);
    }
    ifile.close();

//    return par;
}

} /// namespace jian

