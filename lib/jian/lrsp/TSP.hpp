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

    std::pair<int, int> loop_head_tail(loop *l) {
        int a = t->head->num;
        int b;
        LOOP_EACH(t,
            if (RES->next == NULL) b = RES->num;
        );
        return {a, b};
    }

    void set_pts(loop *l, pts_t &pts) {
        if (l != NULL) {
            auto p = loop_head_tail(l);
            if (p.second - p.first < _max_len) {
                pts.push_back(l);
            } else {
                for (loop *t = l->son; t != NULL; t = t->brother) {
                    set_pts(t, pts);
                }
            }
        }
    }

    Chain *loop_pred(loop *l) {
        if (l != NULL) {
            Chain *chain = find_templ(l);
            auto it = l->hinges.begin();
            for (loop *t = l->son; t != NULL; t = t->brother) {
                auto p = loop_head_tail(t);
                if (p.second - p.first < _max_len) {
                    splice(chain, assemble_mc(t), *it);
                } else {
                    splice(chain, loop_pred(t), *it);
                }
                it++;
            }
            return chain;
        } else {
            return NULL;
        }
    }

    Chain *assemble_mc(loop *l) {
        Chain *chain = assemble(l);
        return mc(chain, l);
    }

    Chain *assemble(loop *l) {
        if (l != NULL) {
            Chain *chain = find_templ(l);
            auto it = l->hinges.begin();
            for (loop *t = l->son; t != NULL; t = t->brother) {
                splice(chain, assemble(t), *it);
                it++;
            }
            return chain;
        } else {
            return NULL;
        }
    }

    Chain *splice(Chain *chain1, Chain *chain2, const hinge_t &hinge) {
        Mat *m1 = mat_chain(chain1, hinge);
        Mat *m2 = mat_chain(chain2);
        delete m1;
        delete m2;
        return chain;
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
        Chain *chain = loop_pred(l);
        chain = mc(chain, l);
        std::cout << *chain << std::endl;
//        traverse(l);
//        pts_t pts;
//        set_pts(l, pts);
//        print_pts(pts);
        free_ss_tree(l);
    }
};

} // namespace lrsp
} // namespace jian

