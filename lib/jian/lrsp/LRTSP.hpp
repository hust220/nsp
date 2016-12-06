#pragma once

#include <thread>
#include <list>
#include <string>
#include "../nuc2d/SSE.hpp"
#include "../nuc2d/SSTree.hpp"
#include "../nuc3d/Assemble.hpp"
#include "../nuc2d/NASS.hpp"
#include "../nuc3d/BuildLoopDG.hpp"
#include "../dhmc/DHMC.hpp"

BEGIN_JN
namespace lrsp {

class LRTSP {
public:
    using pts_t = std::list<SSE *>;
    using record_t = nuc3d::record_t;
    using records_t = nuc3d::records_t;
    using map_records_t = std::map<SSE *, records_t>;
    using map_templates_t = std::map<SSE *, Chain *>;

    map_records_t    _helix_records;
    map_templates_t  _helix_templates;
    map_templates_t  _helix_templates_mc;
    map_records_t    _loop_records;
    map_templates_t  _loop_templates;
    map_templates_t  _loop_templates_mc;

    std::list<Chain *> _gc_chain;
    std::list<std::thread> threads;

    S _name;
    S _seq;
    S _ss;
    S _out;
    S _traj;
    int _max_len = 400;

    LRTSP(const Par &par) {
        par.set(_name, "name", "job");
        par.set(_seq, "seq");
        par.set(_ss, "ss");
        if (par.has("out")) {
            par.set(_out, "out");
        } else {
            std::cout << "Please input the output file name: ";
            std::cin >> _out;
        }
        par.set(_traj, "traj");
    }

    ~LRTSP() {
        for (auto && chain : _gc_chain) {
            delete chain;
        }
    }

    int loop_length(SSE *l) {
        if (l != NULL) {
            auto p = loop_head_tail(l);
            return p.second - p.first + 1;
        } else {
            return 0;
        }
    }

    void set_templates(SSE *l, bool &b_NR) {
        if (l != NULL) {
			auto p = loop_head_tail(l);
            if (l->has_loop() && p.second - p.first < _max_len) {
                set_templates_all(l);
            } else {
                if (l->has_helix()) {
                    set_helix_templates(l, b_NR);
                }
                if (l->has_loop()) {
                    set_loop_templates(l, b_NR);
                }
                for (SSE *t = l->son; t != NULL; t = t->brother) {
                    set_templates(t, b_NR);
                }
            }
        }
    }

    void set_fixed_ranges(SSE *l, fixed_ranges_t &fixed_ranges) {
        if (l != NULL) {
            if (l->has_loop()) {
                auto p = loop_head_tail(l);
                if (p.second - p.first < _max_len) {
                    fixed_ranges.push_back({p.first, p.second});
                } else {
                    for (SSE *t = l->son; t != NULL; t = t->brother) {
                        set_fixed_ranges(t, fixed_ranges);
                    }
                }
            }
        }
    }

    static void lrtsp_chain_refine(Chain *chain, SSE *l, S log) {
        log_file(log);
        chain_refine<CGpsb>(*chain, l, {});
    }

    void set_templates_all(SSE *l) {
        LOGI << "set templates all..." << std::endl;
        bool b_NR = false; // NR: need_refining
        set_templates_mc(l, b_NR);
        print_templates(_loop_templates_mc, "loop templates mc");
        position_templates_mc(NULL, l);
        Chain *chain = make_chain();
        chain->resize(loop_length(l));
        auto p = loop_head_tail(l);
        assemble_templates_mc(*chain, l, p.first);
        if (b_NR) {
            std::ostringstream stream;
            stream << _out << '.' << threads.size() + 1 << ".log";
            threads.push_back(std::thread(lrtsp_chain_refine, chain, l, stream.str()));
//            nuc3d::chain_refine<nuc3d::mc::MCSM>(*chain, l);
        }
        _loop_templates[l] = chain;
        LOGI << l << ' ' << chain->size() << std::endl;
    }

    void set_templates_mc(SSE *l, bool &b_NR) {
        if (l != NULL) {
            if (l->has_helix()) {
                set_helix_templates_mc(l, b_NR);
            }
            if (l->has_loop()) {
                set_loop_templates_mc(l, b_NR);
            }
            for (SSE *t = l->son; t != NULL; t = t->brother) {
                set_templates_mc(t, b_NR);
            }
        }
    }

    Chain *make_chain() {
        Chain *chain = new Chain;
        _gc_chain.push_back(chain);
        return chain;
    }

    void set_helix_templates(SSE *l, bool &b_NR) {
        if (_helix_records[l].empty()) {
            _helix_templates[l] = NULL;
        } else {
            _helix_templates[l] = make_chain();
            nuc3d::chain_read_record(*(_helix_templates[l]), _helix_records[l][0]);
        }
    }

    void set_loop_templates(SSE *l, bool &b_NR) {
        if (_loop_records[l].empty()) {
            b_NR = true;
            _loop_templates[l] = build_chain_dg(l->seq(), NASS::lower_ss(l->ss()));
        } else {
            _loop_templates[l] = make_chain();
            nuc3d::chain_read_record(*(_loop_templates[l]), _loop_records[l][0]);
        }
    }

    void set_helix_templates_mc(SSE *l, bool &b_NR) {
        if (_helix_records[l].empty()) {
            _helix_templates_mc[l] = NULL;
        } else {
            _helix_templates_mc[l] = make_chain();
            nuc3d::chain_read_record(*(_helix_templates_mc[l]), _helix_records[l][0]);
        }
    }

    void set_loop_templates_mc(SSE *l, bool &b_NR) {
        if (_loop_records[l].empty()) {
            b_NR = true;
            _loop_templates_mc[l] = build_chain_dg(l->seq(), NASS::lower_ss(l->ss()));
        } else {
            _loop_templates_mc[l] = make_chain();
            nuc3d::chain_read_record(*(_loop_templates_mc[l]), _loop_records[l][0]);
        }
    }

    void position_templates(Mat *m, SSE *l) {
        if (l != NULL) {
            position_hairpin(m, _helix_templates[l], _loop_templates[l]);
			auto p = loop_head_tail(l);
            if (l->has_loop() && p.second - p.first < _max_len) {
                // pass
            } else if (_loop_templates[l] != NULL) {
                auto it = l->hinges.begin();
                for (SSE *t = l->son; t != NULL; t = t->brother) {
                    m = chain_mat(*(_loop_templates[l]), {it->first-1, it->first, it->second, it->second+1});
                    position_templates(m, t);
                    it++;
                }
            }
        }
    }

    void position_templates_mc(Mat *m, SSE *l) {
        if (l != NULL) {
            position_hairpin(m, _helix_templates_mc[l], _loop_templates_mc[l]);
            if (_loop_templates_mc[l] != NULL) {
                auto it = l->hinges.begin();
                for (SSE *t = l->son; t != NULL; t = t->brother) {
                    m = chain_mat(*(_loop_templates_mc[l]), {it->first-1, it->first, it->second, it->second+1});
                    position_templates_mc(m, t);
                    it++;
                }
            }
        }
    }

    void position_hairpin(Mat *m, Chain *chain_helix, Chain *chain_loop) {
        Chain *chain;
        chain = chain_helix;
        if (chain != NULL) {
            int len = chain->size();
            if (m != NULL) {
                position_template(m, chain);
            }
            m = chain_mat(*chain, {len/2-2, len/2-1, len/2, len/2+1});
        }
        chain = chain_loop;
        if (chain != NULL && m != NULL) {
            position_template(m, chain);
        }
    }

    void position_template(Mat *m, Chain *chain) {
        if (m != NULL) {
            int len = chain->size();
            Mat *mat = chain_mat(*chain, {0, 1, len-2, len-1});
			geom::Superposition<double> sp(*mat, *m);
            for (auto && res : *chain) {
                for (auto && atom : res) {
					sp.apply(atom);
                }
            }
            delete m;
            delete mat;
        }
    }

    Mat *chain_mat(const Chain &chain, const std::list<int> &list) {
        static std::set<std::string> _suppos_atoms{"C5*", "O3*", "C1*"};

        Mat *mat = new Mat(_suppos_atoms.size() * list.size(), 3);
        int index = 0;
        for (int n = 0; n < chain.size(); n++) {
            if (std::find(list.begin(), list.end(), n) != list.end()) {
                for (auto && atom : chain[n]) {
                    if (_suppos_atoms.find(atom.name) != _suppos_atoms.end()) {
                        for (int i = 0; i < 3; i++) {
                            (*mat)(index, i) = atom[i];
                        }
                        index++;
                    }
                }
            }
        }
        return mat;
    }

    void assemble_templates(Chain &chain, SSE *l) {
        LOGI << "assemble loop " << l << std::endl;
        if (l != NULL) {
			auto p = loop_head_tail(l);
            if (l->has_loop() && p.second - p.first < _max_len) {
                // pass
                chain_from_template(chain, l, _loop_templates[l]);
            } else {
                chain_from_template(chain, l, &_helix_templates, &_loop_templates, 0);
                for (SSE *t = l->son; t != NULL; t = t->brother) {
                    assemble_templates(chain, t);
                }
            }
        }
    }

    void assemble_templates_mc(Chain &chain, SSE *l, int beg) {
        LOGI << "assembling templates mc..." << std::endl;
        if (l != NULL) {
            chain_from_template(chain, l, &_helix_templates_mc, &_loop_templates_mc, beg);
            if (l->has_loop()) {
                for (SSE *t = l->son; t != NULL; t = t->brother) {
                    assemble_templates_mc(chain, t, beg);
                }
            }
        }
    }

    void chain_from_template(Chain &chain, SSE *l, Chain *chain_loop) {
        auto p = loop_head_tail(l);
        int n = 0;
        for (int i = p.first; i <= p.second; i++) {
            chain[i] = chain_loop->at(n);
            n++;
        }
    }

    void chain_from_template(Chain &chain, SSE *l, map_templates_t *helix_templates, map_templates_t *loop_templates, int beg) {
        if (l != NULL) {
            std::cout << "beg: " << beg << std::endl;
            l->print();
            if (l->has_loop()) {
                Chain *loop_chain = loop_templates->at(l);
                if (loop_chain == NULL) {
                    throw "chain_from_template error! no loop chain";
                }
                int i = 0;
                for (res *r = l->head; r != NULL; r = r->next) {
                    if (r->type != '&') {
                        chain[r->num - 1 - beg] = loop_chain->at(i);
                        i++;
                    }
                }
                LOGI << std::endl;
            }
            if (l->has_helix()) {
                int len = l->s.len();
                Chain *helix_chain = helix_templates->at(l);
                if (helix_chain == NULL) {
                    throw "chain_from_template error! no helix chain";
                }
                HELIX_EACH(l->s,
                    LOGI << BP->res1.num - 1 - beg << '-' << N_BP << ' ';
                    LOGI << BP->res2.num - 1 - beg << '-' << 2*len-1-N_BP << std::endl;
                    chain[BP->res1.num - 1 - beg] = (*helix_chain)[N_BP];
                    chain[BP->res2.num - 1 - beg] = (*helix_chain)[2*len-1-N_BP];
                );
            }
        }
    }

    void find_records(SSE *l) {
        if (l != NULL) {
            l->print();
            nuc3d::find_loop_records(l, _loop_records[l]);
            nuc3d::find_helix_records(l, _helix_records[l]);
            for (SSE *t = l->son; t != NULL; t = t->brother) {
                find_records(t);
            }
        }
    }

    void print_templates(const map_templates_t &m, S name) {
        LOGI << name << std::endl;
        for (auto && pair : m) {
            if (pair.second != NULL) {
                LOGI << pair.first << ' ' << (pair.second) << ' ' << pair.second->size() << std::endl;
            }
        }
    }

    void print_templates() {
        print_templates(_helix_templates, "helix templates");
        print_templates(_loop_templates, "loop templates");
        print_templates(_helix_templates_mc, "helix templates of mc");
        print_templates(_loop_templates_mc, "loop templates of mc");
    }

    void print_records() {
        LOGI << "helix records" << std::endl;
        for (auto && pair : _helix_records) {
            LOGI << pair.first << ' ' << pair.second.size() << std::endl;
        }
        LOGI << "loop records" << std::endl;
        for (auto && pair : _loop_records) {
            LOGI << pair.first << ' ' << pair.second.size() << std::endl;
        }
    }

    void init_templates(SSE *l) {
		for (auto &&sse : pretree(l)) {
			SSE *_l = &sse;
			_helix_templates[_l] = NULL;
			_loop_templates[_l] = NULL;
			_helix_templates_mc[_l] = NULL;
			_loop_templates_mc[_l] = NULL;
		}
    }

    void pred() {
        LOGI << "# constructing ss tree..." << std::endl;
        SSE *l = ss_tree(_seq, _ss);
        LOGI << "# initializing templates..." << std::endl;
        init_templates(l);
        LOGI << "# searching for records..." << std::endl;
        find_records(l);
        LOGI << "# printing records..." << std::endl;
        print_records();
        bool b_NR;
        LOGI << "# setting templates..." << std::endl;
        set_templates(l, b_NR);
        LOGI << "# waiting for refining..." << std::endl;
        for (auto && thread : threads) {
            thread.join();
        }
        LOGI << "# printing templates..." << std::endl;
        print_templates();
        LOGI << "# positioning templates..." << std::endl;
        position_templates(NULL, l);
        LOGI << "# assembling templates..." << std::endl;
        Chain *chain = make_chain();
        chain->resize(_seq.size());
        assemble_templates(*chain, l);
        if (b_NR) {
            LOGI << "# refining..." << std::endl;
            fixed_ranges_t fixed_ranges;
            set_fixed_ranges(l, fixed_ranges);
            chain_refine<CGpsb>(*chain, l, fixed_ranges, _traj);
        }
        mol_write(*chain, _out);
        LOGI << "# free ss tree..." << std::endl;
        free_ss_tree(l);
    }
};

} // namespace lrsp
END_JN

