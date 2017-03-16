#pragma once

#include <jian/utils/traits.hpp>
#include "../nuc2d.hpp"

BEGIN_JN

class TemplRec {
public:
    Str name;
    Str seq;
    Str ss;
    Str family;
    Str src;
    Int len = 0;
    Int type;
    Num score = 0;
};

Bool set_loop_rec(TemplRec &, const Str &);

Bool set_helix_rec(TemplRec &, const Str &);

Deque<TemplRec> find_helix_templates(Str seq, Str src = "");

Deque<TemplRec> find_loop_templates(Str seq, Str ss, Str src = "");

END_JN

