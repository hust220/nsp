#include "../utils/string.hpp"
#include "TemplRec.hpp"

namespace jian {

void templ_rec_set_src(TemplRec &rec) {
    int pos = rec._name.find_first_of("-");
    if (pos == str_t::npos) {
        rec._src = jian::upper(rec._name);
    } else {
        rec._src = jian::upper(rec._name.substr(0, pos));
    }
}

bool set_loop_rec(TemplRec &rec, const std::string &s) {
    tokenize_v v;
    jian::tokenize(s, v, " ");
    if (v.size() != 5) {
        return false;
    } else {
        rec._name = v[0];
        templ_rec_set_src(rec);
        rec._type = JN_INT(v[1]);
        rec._seq = v[2];
        rec._ss = v[3];
        rec._family = v[4];
        return true;
    }
}

bool set_helix_rec(TemplRec &rec, const std::string &s) {
    tokenize_v v;
    jian::tokenize(s, v, " ");
    if (v.size() != 5) {
        return false;
    } else {
        rec._name = v[0];
        templ_rec_set_src(rec);
        rec._len = JN_INT(v[1]);
        rec._seq = v[2];
        rec._ss = v[3];
        rec._family = v[4];
        return true;
    }
}

}

