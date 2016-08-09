#pragma once

#include <list>
#include <string>
#include "../nuc2d/loop.hpp"
#include "../nuc2d/SSTree.hpp"

namespace jian {
namespace lrsp {

class TSP {
public:
    using pts_t = std::list<loop *>;

    std::string _seq;
    std::string _ss;
    int _max_len = 400;

    TSP(const Par &par) {
        par.set(_seq, "seq");
        par.set(_ss, "ss");
    }

    void print_pts(const pts_t & pts) {
        for (auto && pt : pts) {
            LOGI << pt << ' ';
        }
        LOGI << std::endl;
    }

    void set_pts(loop *l, pts_t &pts) {
        if (l != NULL) {
            int a = l->head->num;
            int b;
            LOOP_EACH(l,
                if (RES->next == NULL) b = RES->num;
            );
            if (b - a < _max_len) {
                pts.push_back(l);
            } else {
                for (loop *t = l->son; t != NULL; t = t->brother) {
                    set_pts(t, pts);
                }
            }
        }
    }

    void traverse(loop *l) {
        if (l != NULL) {
            l->print();
            int a = l->head->num;
            int b;
            LOOP_EACH(l,
                if (RES->next == NULL) b = RES->num;
            );
            if (b - a < _max_len) {
                LOGI << "pass: " << l << std::endl;
            } else {
                for (loop *t = l->son; t != NULL; t = t->brother) {
                    traverse(t);
                }
            }
        }
    }

    void pred() {
        loop *l = ss_tree(_seq, _ss);
        traverse(l);
        pts_t pts;
        set_pts(l, pts);
        print_pts(pts);
        free_ss_tree(l);
    }
};

} // namespace lrsp
} // namespace jian

