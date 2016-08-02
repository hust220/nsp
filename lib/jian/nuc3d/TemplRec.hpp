#pragma once

#include <string>

namespace jian {

class TemplRec {
public:
    std::string _name;
    std::string _seq;
    std::string _ss;
    std::string _family;
    std::string _src;
    int _len = 0;
    int _type;
    double _score = 0;
};

bool set_loop_rec(TemplRec &, const std::string &);
bool set_helix_rec(TemplRec &, const std::string &);

}

