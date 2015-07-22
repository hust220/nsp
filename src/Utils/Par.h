#ifndef PAR_H
#define PAR_H

#include "std.h"
#include "MLib.h"

namespace jian {

class Par: public std::map<std::string, std::vector<std::string>> {
public:
    Par() {}
    Par(int argc, char **argv) {
        read(argc, argv);
//        _par = read(argc, argv);
    }
    Par(string str) {
        read(str);
//        _par = read(str);
    }
    void read(int, char **);
    void read(string);
//    static map<string, vector<string>> read(int, char **);
//    static map<string, vector<string>> read(string);

//    map<string, vector<string>> _par;
};

} /// namespace jian

#endif
