#pragma once

#include <string>

BEGIN_JN

class TemplRec {
public:
    S _name;
    S _seq;
    S _ss;
    S _family;
    S _src;
    int _len = 0;
    int _type;
    double _score = 0;
};

bool set_loop_rec(TemplRec &, const S &);
bool set_helix_rec(TemplRec &, const S &);

}

