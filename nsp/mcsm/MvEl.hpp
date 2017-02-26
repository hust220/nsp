#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../nuc3d/TSP.hpp"
#include "../nuc3d/BuildHelix.hpp"
#include "../nuc3d/transform.hpp"
#include "../nuc3d/TemplRec.hpp"
#include "../cg.hpp"
#include "../cg/ResFrags.hpp"
#include "jian/utils/Env.hpp"
#include "jian/utils/file.hpp"
#include "../mc.hpp"
#include "../pdb.hpp"
#include "jian/geom.hpp"
#include "../nuc2d.hpp"
#include "jian/pp.hpp"

#pragma once

BEGIN_JN

using Frag = Array<int, 2>;

using Frags = Deque<Frag>;

class MvEl;

using MvEls = Deque<MvEl>;

// MvEl: Moving Element
class MvEl {
    public:
        enum Type {
            MVEL_HL, // helix
            MVEL_HP, // hairpin
            MVEL_IL, // internal loop
            MVEL_FG // fragment
        };

        Type type;
        Frags range;

        MvEl(int a, int b, Type t);

        MvEl(int a, int b, int c, int d, Type t);

        MvEl(const Helix &h);

        MvEl(SSTree::El *l, Type t);

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

