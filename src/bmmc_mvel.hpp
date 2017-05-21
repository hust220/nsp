#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "rtsp.hpp"
#include "rtsp_templ_rec.hpp"
#include "cg.hpp"
#include "cg_res_frags.hpp"
#include "env.hpp"
#include "file.hpp"
#include "mc.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "rss.hpp"
#include "pp.hpp"

#pragma once

BEGIN_JN

using Frag = Array<int, 2>;

using Frags = Deque<Frag>;

Frag frag_read(Str s);

Frags frags_read(Str s);

class MvEl;

using MvEls = Deque<MvEl>;

enum {
    MVEL_HL,
    MVEL_IL,
    MVEL_HP,
    MVEL_FG,
    MVEL_ML
};

// MvEl: Moving Element
class MvEl {
    public:
        using Type = Int;

        Type type;
        Frags range;

        MvEl(Type t);

        MvEl(int a, int b, Type t);

        MvEl(int a, int b, int c, int d, Type t);

        MvEl(const Helix &h);

        MvEl(SSTree::El *l, Type t);

        MvEl &add_frag(int a, int b);

        int min() const;

        int max() const;

        bool operator ==(const MvEl &el) const;

        bool operator !=(const MvEl &el) const;

        MvEl *operator +(const MvEl &el) const;

        friend std::ostream &operator <<(std::ostream &, const MvEl &el);

        bool contains(const MvEl &el) const;

        bool nips(const MvEl &el) const;

        bool has(int n) const {
            return std::find_if(range.begin(), range.end(), [&n](const Frag &frag) {
                    return frag[0] <= n && n <= frag[1];
                    }) != range.end();
        }

        bool minmax_has(int n) const {
            int min = 99999;
            int max = -1;
            for (auto && frag : range) {
                if (min > frag[0]) min = frag[0];
                if (max < frag[1]) max = frag[1];
            }
            return min <= n && n <= max;
        }

        static void merge(Deque<MvEl *> &dq);
};

END_JN

