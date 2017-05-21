#include "bmmc_mvel.hpp"

BEGIN_JN

Frag frag_read(Str s) {
    tokenize_v v;
    int beg, end;

    tokenize(s, v, "-");
    beg = JN_INT(v[0]) - 1;
    if (v.size() == 1) {
        end = beg;
    }
    else if (v.size() == 2) {
        end = JN_INT(v[1]) - 1;
    }
    return {beg, end};
}

Frags frags_read(Str str) {
    Frags frags;
    tokenize_v w;

    tokenize(str, w, "+");
    for (auto && s : w) {
        frags.push_back(frag_read(s));
    }
    return frags;
}


MvEl::MvEl(MvEl::Type t) : type(t) {}

MvEl &MvEl::add_frag(int a, int b) {
    range.push_back({a, b});
}

MvEl::MvEl(int a, int b, MvEl::Type t) : type(t) {
    range.push_back({ a, b });
}

MvEl::MvEl(int a, int b, int c, int d, MvEl::Type t) : type(t) {
    range.push_back({ a, b });
    range.push_back({ c, d });
}

MvEl::MvEl(const Helix &h) : type(MVEL_HL) {
    int a, b, c, d;

    a = h.front().res1.num - 1;
    d = h.front().res2.num - 1;
    for (auto && bp : h) {
        if (bp.next == NULL) {
            b = bp.res1.num - 1;
            c = bp.res2.num - 1;
        }
    }
    range.push_back({ a, b });
    range.push_back({ c, d });
}

MvEl::MvEl(SSTree::El *l, MvEl::Type t) : type(t) {
    int a, b, c, d;

    if (t == MVEL_HP) {
        a = l->data.helix.front().res1.num - 1;
        b = l->data.helix.front().res2.num - 1;
        range.push_back({ a, b });
    }
    else if (t == MVEL_IL) {
        a = l->data.helix.front().res1.num - 1;
        b = l->son->data.helix.front().res1.num - 2;
        c = l->son->data.helix.front().res2.num;
        d = l->data.helix.front().res2.num - 1;
        range.push_back({ a, b });
        range.push_back({ c, d });
    }
    else {
        throw "jian::MvEl error!";
    }
}

bool MvEl::operator ==(const MvEl &el) const {
    return type == el.type && range == el.range;
}

bool MvEl::operator !=(const MvEl &el) const {
    return !(*this == el);
}

MvEl *MvEl::operator +(const MvEl &el) const {
    if (el.range.size() == 1) {
        return new MvEl(range[0][0], range[1][1], el.type);
    }
    else if (el.range.size() == 2) {
        return new MvEl(range[0][0], el.range[0][1], el.range[1][0], range[1][1], type);
    }
    else {
        throw "jian::MvEl *jian::MvEl::operator +(const jian::MvEl &el) error!";
    }
}

std::ostream &operator <<(std::ostream &stream, const MvEl &el) {
    stream <<
        (el.type == MVEL_HL ? "Helix" :
         (el.type == MVEL_HP ? "Hairpin" :
          (el.type == MVEL_IL ? "Internal Loop" :
           (el.type == MVEL_FG ? "Fragment" : "Others")))) << ' ';
    for (auto && frag : el.range) {
        stream << '(' << frag[0] << '-' << frag[1] << ')';
    }
    return stream;
}

int MvEl::min() const {
    return std::min_element(range.begin(), range.end(), [](const Frag &f1, const Frag &f2) {
            return f1[0] <= f2[0];
            })->at(0);
}

int MvEl::max() const {
    return std::max_element(range.begin(), range.end(), [](const Frag &f1, const Frag &f2) {
            return f1[1] <= f2[1];
            })->at(1);
}

bool MvEl::contains(const MvEl &el) const {
    return std::all_of(el.range.begin(), el.range.end(), [this](const Frag &f1) {
            return std::any_of(this->range.begin(), this->range.end(), [&f1](const Frag &f) {
                return f[0] <= f1[0] && f[1] >= f1[1];
                });
            });
}

bool MvEl::nips(const MvEl &el) const {
    return range.size() == 2 && range[0][1] + 1 == el.min() && range[1][0] - 1 == el.max();
}

void MvEl::merge(Deque<MvEl *> &dq) {
    //log << "# Merge ranges..." << std::endl;
    Deque<MvEl *> els;
    int flag = 1;
    Map<MvEl *, Bool> m;

    if (dq.empty()) return;

    while (flag != 0) {
        flag = 0;
        m.clear();
        for (auto && el : dq) m[el] = true;

        auto it = dq.begin();
        for (auto it1 = it; it1 != dq.end(); it1++) {
            if (!m[*it1]) continue;
            for (auto it2 = it1 + 1; it2 != dq.end(); it2++) {
                if (!m[*it2]) continue;
                if ((*it1)->type != MVEL_FG && (*it2)->type != MVEL_FG) {
                    MvEl &el1 = *(*it1);
                    MvEl &el2 = *(*it2);
                    if (el1.contains(el2)) {
                        m[*it2] = false;
                        flag++;
                    }
                    else if (el2.contains(el1)) {
                        m[*it1] = false;
                        flag++;
                    }
                    else if (el1.nips(el2)) {
                        m[*it1] = false;
                        m[*it2] = false;
                        els.push_back(el1 + el2);
                        flag++;
                    }
                    else if (el2.nips(el1)) {
                        m[*it1] = false;
                        m[*it2] = false;
                        els.push_back(el2 + el1);
                        flag++;
                    }
                }
            }
        }

        for (auto && el : dq) {
            if (m[el]) {
                els.push_back(el);
            }
            else {
                delete el;
            }
        }
        dq = els;
        els.clear();
    }
}



END_JN
