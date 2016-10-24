#include <array>
#include <algorithm>
#include <mutex>
#include <list>
#include <vector>
#include <map>
#include <string>
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "../utils/log.hpp"

#include "ss_pred.hpp"

namespace jian {

namespace ss_pred_detail {

std::mutex g_mt;

using val_t = double;

struct loop_t {
    int beg, end;
};

using path_t = std::list<loop_t *>;

using inf_t = struct {val_t score; path_t path;};

void print_path(const path_t &path) {
    for (auto && l : path) LOGV << l->beg << '-' << l->end << ' '; LOGV << std::endl;
}

class Pred2D {
public:
    //using sons_t = std::list<son_t *>;
    using mat_t = std::array<std::array<val_t, 4>, 4>;
    using stack_t = std::array<std::array<mat_t, 4>, 4>;

    std::string _lib;

    std::vector<int> _indices;
    std::map<int, std::map<int, inf_t*>> _inf;
    std::list<inf_t *> gc_inf;
    std::list<loop_t *> gc_loop;

    std::array<val_t, 30> _en_hairpin;
    std::array<val_t, 30> _en_internal;
    std::array<val_t, 30> _en_bulge;
    stack_t _en_stack;
    val_t _offset,  _free_base_penalty,  _helix_penalty;

    Pred2D() {
        _lib = Env::lib() + "/RNA/pars/nuc2d/mfold";
        read_loop_dg();
        read_stack_dg();
        read_miscloop_dg();
    }

    template<typename T>
    void print_arr(T && arr) {
        for (auto && i : arr) {
            LOGV << i << std::endl;
        }
    }

    void print_pars() {
        LOGV << "Internal loop parameters: " << std::endl;
        print_arr(_en_internal);
        LOGV << "Bulge loop parameters: " << std::endl;
        print_arr(_en_bulge);
        LOGV << "Hairpin loop parameters: " << std::endl;
        print_arr(_en_hairpin);
        LOGV << "Stacking parameters: " << std::endl;
        for (auto && i : _en_stack) {
            for (auto && j : i) {
                for (auto && k : j) {
                    for (auto && l : k) {
                        LOGV << l << ' ';
                    }
                    LOGV << std::endl;
                }
                LOGV << std::endl;
            }
        }
        LOGV << "Junction parameters: " << std::endl;
        LOGV << _offset << ' ' << _free_base_penalty << ' ' << _helix_penalty << std::endl;
    }

    void read_loop_dg() {
        LOGV << "### Read loop.dg" << std::endl;
        int n; std::vector<std::string> v;
        EACH_LINE((_lib + "/loop.dg").c_str(),
            if (N+1 >= 5) {
                jian::trim(L);
                tokenize(L, v, " ");
                if (v.size() == 4) {
                    n = std::stoi(v[0]) - 1; 
                    _en_internal[n] = (v[1] == "." ? 0 : std::stof(v[1]));
                    _en_bulge[n] = (v[2] == "." ? 0 : std::stof(v[2]));
                    _en_hairpin[n] = (v[3] == "." ? 0 : std::stof(v[3]));
                }
            }
        );
    }

    void read_stack_dg() {
        LOGV << "### Read stack.dg" << std::endl;
        std::vector<std::string> v;
        auto foo = [&](auto &line, int m, int i){
            jian::trim(line); tokenize(line, v, " ");
            if (v.size() == 16) {
                for (int n = 0; n < 4; n++) for (int j = 0; j < 4; j++) {
                    _en_stack[m][n][i][j] = (v[n * 4 + j] == "." ? 0 : std::stod(v[n * 4 + j]));
                }
            }
        };
        EACH_LINE((_lib + "/stack.dg").c_str(),
            int n = N+1;
            if (27 <= n && n <=30) foo(L, 0, n - 27);
            else if (41 <= n && n <= 44) foo(L, 1, n - 41);
            else if (56 <= n && n <= 59) foo(L, 2, n - 56);
            else if (70 <= n && n <= 73) foo(L, 3, n - 70);
        );
    }

    void read_miscloop_dg() {
        LOGV << "### Read miscloop.dg" << std::endl;
        int n = 0;
        std::vector<std::string> v;
        EACH_LINE((_lib + "/miscloop.dg").c_str(),
            if (std::regex_search(L, std::regex("offset,  free base penalty,  helix penalty"))) {
                n++;
            } else if (n == 1) {
                jian::trim(L);
                tokenize(L, v, " ");
                _offset = std::stod(v[0]);
                _free_base_penalty = std::stod(v[1]);
                _helix_penalty = std::stod(v[2]);
                n++;
            }
        );
    }

    void init() {
        _inf.clear();
        for (auto && l : gc_loop) delete l; gc_loop.clear();
        for (auto && l : gc_inf) delete l; gc_inf.clear();
    }

    ~Pred2D() {
        init();
    }

    void seq_to_indices(const std::string &seq) {
        std::map<char, int> temp_map{{'A', 0}, {'C', 1}, {'G', 2}, {'U', 3}, {'X', -1}};
        _indices.resize(seq.size());
        std::transform(seq.begin(), seq.end(), _indices.begin(), [&](char c){return temp_map[c];});
        for (auto && i : _indices) LOGV << i << ' '; LOGV << std::endl;
    }

    std::string operator ()(const std::string &seq) {
        std::lock_guard<std::mutex> gd(g_mt);
        print_pars();
        init();
        seq_to_indices(seq);
        std::string ss(seq.size(), '.');
        loop_inf(-1, seq.size());
        set_loop(ss, -1, seq.size());
        return ss;
    }

    void set_loop(std::string &s, int beg, int end) {
        if (_inf.count(beg) && _inf[beg].count(end)) {
            if (!(_inf[beg][end]->path.empty())) {
                for (auto && loop : _inf[beg][end]->path) {
                    s[loop->beg] = '(';
                    s[loop->end] = ')';
                    set_loop(s, loop->beg, loop->end);
                }
            }
        }
    }

    val_t loop_score(int beg, int end, const path_t &s) {
        int l_d = s.size();
        int l_s = end - beg - 1;
        for (auto &&l : s) l_s -= (l->end - l->beg + 1);
        val_t score = 0;
        if (l_d == 0) {
            if (l_s > 0) score += en_hairpin(l_s);
        } else if (l_d == 1) {
            if (l_s == 0) {
                if (beg != -1 && end != _indices.size()) {
                    score += en_stack(beg, end, s.front()->beg, s.front()->end);
                }
            } else if (s.front()->beg - beg == 1 || end - s.front()->end == 1) score += en_bulge(l_s);
            else score += en_internal(l_s);
        } else {
            score += en_junction(l_s, l_d);
        }
        for (auto && l : s) {
            inf_t *inf = loop_inf(l->beg, l->end);
            score += (inf == NULL ? 0 : inf->score);
            //score += loop_inf(l->beg, l->end)->score;
        }
        return score;
    }

    val_t en_hairpin(int l_s) {
        if (l_s <= 30) {
            return _en_hairpin[l_s-1];
        } else {
            return _en_hairpin[29]+0.05*(l_s-30);
        }
    }

    val_t en_stack(int i, int j, int k, int l) {
//        auto check = [](int n){
//            if (n < 0 || n > 3) return false;
//            return true;
//        };
        int a = _indices[i];
        int b = _indices[j];
        int c = _indices[k];
        int d = _indices[l];
//        if ( check(a) && check(b) && check(c) && check(d)) {
//            
//        } else {
//            LOGI << i << ' ' << j << ' ' << k << ' ' << l << std::endl;
//            LOGI << a << ' ' << b << ' ' << c << ' ' << d << std::endl;
//        }
        if (a == -1 || b == -1 || c == -1 || d == -1) return 0;
        else return 2 * _en_stack[a][b][c][d];
    }

    val_t en_bulge(int l_s) {
        if (l_s <= 30) {
            return _en_bulge[l_s-1];
        } else {
            return _en_bulge[29]+0.05*(l_s-30);
        }
    }

    val_t en_internal(int l_s) {
        if (l_s <= 30) {
            return _en_internal[l_s-1];
        } else {
            return _en_internal[29]+0.05*(l_s-30);
        }
    }

    val_t en_junction(int l_s, int l_d) {
        return 0.2 * (_offset +  _free_base_penalty * l_s +  _helix_penalty * l_d);
    }

    inf_t *loop_inf(int beg, int end) {
        if (end - beg < 4) {
            return NULL;
        } else if (_inf.count(beg) && _inf[beg].count(end)) {
            return _inf[beg][end];
        } else {

            inf_t *inf;
            inf_t *inf1 = helper_loop_inf(beg+1, end);
            inf_t *inf2 = loop_inf(beg+1, end);
            if (inf2 == NULL) {
                inf = inf1;
            } else {
                double score = loop_score(beg, end, inf2->path);
                if (inf1 != NULL && inf1->score < score) {
                    inf = inf1;
                } else {
                    inf = new inf_t {score, inf2->path};
                    gc_inf.push_back(inf);
                }
            }
            _inf[beg][end] = inf;
            return inf;
        }
    }

    inf_t *helper_loop_inf(int i, int end) {
        LOGV << "### Find sons of left-fixed loop (" << i << ", " << end-1 << ")" << std::endl;

        inf_t *best_inf = new inf_t {99999, path_t{}};
        gc_inf.push_back(best_inf);
        for (int j = i + 4; j < end; j++) {
            if (is_pair(_indices[i], _indices[j])) {
                loop_t *l = new loop_t{i, j};
                gc_loop.push_back(l);
                inf_t *inf = loop_inf(j, end);
                path_t path {l};
                if (inf != NULL) {
                    for (auto && n : inf->path) {
                        path.push_back(n);
                    }
                }
                double score = loop_score(i-1, end, path);
                if (score < best_inf->score) {
                    best_inf->score = score;
                    best_inf->path = path;
                }
            }
        }
        return best_inf;
    }

    bool is_pair(int i, int j) {
        static std::map<int, int> m{{0, 1}, {1, 2}, {2, 4}, {3, 8}, {-1, 16}};
        auto n = m[i] | m[j];
        return n == 9 || n == 6 || n == 12;
    }

};

Pred2D pred;

} // namespace_ss_pred_detail

std::string ss_pred(const std::string &seq) {
    return ss_pred_detail::pred(seq);
}

} // namespace jian

