#pragma once

namespace jian {
namespace nuc2d {

class Pred2D {
public:
    using val_t = double;
    using Loop = struct {int beg; int end;};
    using Sons = std::list<Loop>;
    using Inf = struct {val_t score; Sons sons;};
    using Mat = std::array<std::array<val_t, 4>, 4>;
    using MatStack = std::array<std::array<Mat, 4>, 4>;

    std::string _lib;

    std::vector<int> _indices;
    std::map<int, std::map<int, Inf>> _inf;
    std::map<int, std::map<int, std::list<Sons>>> _cache;
    std::map<int, std::map<int, std::list<Sons>>> _cache2;

    std::array<val_t, 30> _en_hairpin;
    std::array<val_t, 30> _en_internal;
    std::array<val_t, 30> _en_bulge;
    MatStack _en_stack;
    val_t _offset,  _free_base_penalty,  _helix_penalty;

    Pred2D(const std::string &lib = "/home/jian/share/mfold") : _lib(lib) {
        log("# Predict 2D structure\n", 
            "## Read Parameters\n");
        read_loop_dg();
        read_stack_dg();
        read_miscloop_dg();
    }

    void read_loop_dg() {
        log("### Read loop.dg\n");
        int n; std::vector<std::string> v;
        file::each_line(_lib + "/loop.dg", [&](auto &line, int num_line){
            if (num_line >= 5) {
                boost::trim(line); tokenize(line, v, " ");
                if (v.size() == 4) {
                    n = std::stoi(v[0]) - 1; 
                    _en_internal[n] = (v[1] == "." ? 0 : std::stof(v[1]));
                    _en_bulge[n] = (v[2] == "." ? 0 : std::stof(v[2]));
                    _en_hairpin[n] = (v[3] == "." ? 0 : std::stof(v[3]));
                }
            }
            return true;
        });
    }

    void read_stack_dg() {
        log("### Read stack.dg\n");
        std::vector<std::string> v;
        auto foo = [&](auto &line, int m, int i){
            boost::trim(line); tokenize(line, v, " ");
            if (v.size() == 16) {
                for (int n = 0; n < 4; n++) for (int j = 0; j < 4; j++) {
                    _en_stack[m][n][i][j] = (v[n * 4 + j] == "." ? 0 : std::stod(v[n * 4 + j]));
                }
            }
        };
        file::each_line(_lib + "/stack.dg", [&](auto &line, int num_line){
            if (27 <= num_line && num_line <=30) foo(line, 0, num_line - 27);
            else if (41 <= num_line && num_line <= 44) foo(line, 1, num_line - 41);
            else if (56 <= num_line && num_line <= 59) foo(line, 2, num_line - 56);
            else if (70 <= num_line && num_line <= 73) foo(line, 3, num_line - 70);
            return true;
        });
    }

    void read_miscloop_dg() {
        log("### Read miscloop.dg\n");
        int n = 0; std::vector<std::string> v;
        file::each_line(_lib + "/miscloop.dg", [&](auto &line, int num_line){
            if (std::regex_search(line, std::regex("offset,  free base penalty,  helix penalty"))) n++;
            else if (n == 1) {
                boost::trim(line); tokenize(line, v, " ");
                _offset = std::stod(v[0]); _free_base_penalty = std::stod(v[1]); _helix_penalty = std::stod(v[2]);
                return false;
            }
            return true;
        });
    }

    val_t operator ()(const std::string &seq) {
        log("## Initialize\n");
        _inf.clear(); _cache.clear(); _cache2.clear();
        log("## Convert sequence to indices\n");
        std::map<char, int> temp_map{{'A', 0}, {'C', 1}, {'G', 2}, {'U', 3}};
        std::transform(seq.begin(), seq.end(), std::back_inserter(_indices), [&](char c){return temp_map[c];});
        for (auto && i : _indices) std::cout << i << ' '; std::cout << std::endl;
        score_loop(-1, seq.size());
        std::string ss(seq.size(), '.');
        set_loop(ss, -1, seq.size());
        std::cout << ss << std::endl;
        return _inf[-1][seq.size()].score;
    }

    void set_loop(std::string &s, int beg, int end) {
        if (beg != -1) {
            s[beg] = '('; s[end] = ')';
        }
        if (_inf.count(beg) && _inf[beg].count(end)) {
            for (auto && loop : _inf[beg][end].sons) {
                s[loop.beg] = '('; s[loop.end] = ')';
                set_loop(s, loop.beg, loop.end);
            }
        }
    }

    val_t score_loop(int beg, int end) {
        if (_inf.count(beg) && _inf[beg].count(end)) return _inf[beg][end].score;
        log("## Get the score of loop(", beg, ", ", end, ")\n");
        auto all_sons = find_sons(beg, end); 
        if (all_sons.empty()) {
            val_t score = helper_score_loop(beg, end, Sons());
            _inf[beg][end] = {score, Sons()};
            log("The score of loop(", beg, ", ", end, "): ", score, "\n");
            return helper_score_loop(beg, end, Sons());    
        } else {
            std::vector<val_t> scores(all_sons.size());
            std::transform(all_sons.begin(), all_sons.end(), scores.begin(), [&](const auto &s){return this->helper_score_loop(beg, end, s);});
            auto it = std::min_element(scores.begin(), scores.end());
            _inf[beg][end] = {*it, *std::next(all_sons.begin(), std::distance(scores.begin(), it))};
            log("The score of loop(", beg, ", ", end, "): ", *it, "\n");
            return *it;
        }
    }

    void print_sons(const Sons &s) {
        for (auto && loop : s) std::cout << loop.beg << '-' << loop.end << ' '; std::cout << "\n";
    }

    val_t helper_score_loop(int beg, int end, const Sons &s) {
        //log("### Helper of get the score of loop(", beg, ", ", end, ")\n");
        int l_s = end - beg - 1; for (auto &&l : s) l_s -= (l.end - l.beg + 1);
        val_t score = 0; int l_d = s.size();
        if (l_d == 0) {
            if (l_s > 0) score += en_hairpin(l_s);
        } else if (l_d == 1) {
            if (l_s == 0 && beg != -1) score += _en_stack[_indices[beg]][_indices[end]][_indices[s.front().beg]][_indices[s.front().end]];
            else if (s.front().beg - beg == 1 || end - s.front().end == 1) score += en_bulge(l_s);
            else score += en_internal(l_s);
        } else score += en_junction(l_s, l_d);
        score += std::accumulate(s.begin(), s.end(), 0, [&](int sum, const auto &l){return sum + this->score_loop(l.beg, l.end);});
        return score;
    }

    val_t en_hairpin(int l_s) {
        return _en_hairpin[l_s];
    }

    val_t en_stack() {}

    val_t en_bulge(int l_s) {
        return _en_bulge[l_s];
    }

    val_t en_internal(int l_s) {
        return _en_internal[l_s];
    }

    val_t en_junction(int l_s, int l_d) {
        return _offset +  _free_base_penalty * l_s +  _helix_penalty * l_d;
    }

    std::list<Sons> find_sons(int beg, int end) {
        //log("### Find the sons of the loop(", beg, ", ", end, ")\n");
        if (_cache.count(beg) && _cache[beg].count(end)) return _cache[beg][end];
        std::list<Sons> ls;
        beg++; end--;
        for (int i = beg; i < end; i++) {
            ls.splice(ls.end(), helper_find_sons(i, end));
        }
        //for (auto && sons : ls) print_sons(sons);
        _cache[beg-1][end+1] = ls;
        return ls;
    }

    std::list<Sons> helper_find_sons(int beg, int end) {
        //log("#### Helper of find the sons of the loop(", beg, ", ", end, ")\n");
        if (_cache2.count(beg) && _cache2[beg].count(end)) return _cache2[beg][end];
        std::list<Sons> ls;
        for (int i = beg + 4; i <= end; i++) {
            if (is_pair(_indices[beg], _indices[i])) {
                auto sons = find_sons(i + 1, end);
                if (sons.empty()) ls.push_back({{beg, i}});
                else for (auto &&s : sons) {
                    Sons ls2{{beg, i}}; ls2.splice(ls2.end(), s); ls.push_back(ls2);
                }
            }
        }
        _cache2[beg][end] = ls;
        return ls;
    }

    bool is_pair(int i, int j) {
        static std::map<int, unsigned int> m{{0, 1}, {1, 2}, {2, 4}, {3, 8}};
        auto n = m[i] | m[j]; return n == 9 || n == 6 || n == 12;
    }

};

} // namespace nuc2d
} // namespace jian

