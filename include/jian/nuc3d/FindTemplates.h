#ifndef JIAN_NUC3D_FINDTEMPLATES_H
#define JIAN_NUC3D_FINDTEMPLATES_H

#include "../util/std.h"
#include "../nuc2d/util.h"
#include "../nuc2d/N2D.h"
#include "TemplRec.h"
#include "JobInf.h"

namespace jian {
namespace nuc3d {

class FindTemplates : public virtual JobInf {
public:
    typedef std::tuple<nuc2d::loop *, int, int, int> Row;

    std::map<nuc2d::loop *, std::pair<std::deque<TemplRec>, std::deque<TemplRec>>> _records;
    std::map<nuc2d::loop *, std::pair<Model, Model>> _templates;
    std::list<Row> _loop_nums_table;

    FindTemplates() {}

    void find_templates() {
        find_records(_n2d.head);
        print_records();
    }

    void print_records() {
        log("Records searching results:\n");
        for (auto &&Row: _loop_nums_table) {
            nuc2d::loop *l;
            int type, num_loops, num_helices;
            std::tie(l, type, num_loops, num_helices) = Row;
            log("loop(", l, "):\n", "Helix: ", l->s.seq(), ' ', l->s.ss(), ' ', num_helices, '\n');
            log("Loop: ", l->seq(), ' ', l->ss(), ' ', num_loops, "\n\n");
        }
    }

    void find_records(nuc2d::loop *l) {
        if (l == NULL) return;
        find_loop_records(l);
        find_helix_records(l);
        find_records(l->son);
        find_records(l->brother);
    }

    void find_loop_records(nuc2d::loop *l) {
        if (l->empty()) {
            _loop_nums_table.push_back(std::make_tuple(l, 0, 0, 0));    
            return;
        }

        std::string seq = l->seq(), ss = l->ss(), p_ss = nuc2d::pure_ss(ss), lower_ss = nuc2d::lower_ss(p_ss, 1), family = _family;
        int num_sons = l->num_sons();

        std::string info_file = _lib + "/" + _type + "/" + "records/" + (l->is_open() ? "open_" : "") + "loop";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::nuc3d::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        int num = 0;
        while (ifile >> templ_rec._name >> templ_rec._type >> templ_rec._seq >> templ_rec._ss >> templ_rec._family) {
            if (_strategy == "loose" and templ_rec._type == num_sons and num_sons >= (l->is_open() ? 0 : 2)) {
                templ_rec._score = 0;
                if (templ_rec._name.substr(0, 4) == _name.substr(0, 4)) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                _records[l].first.push_back(templ_rec);
                num++;
                if (not _source_pdb.empty() and templ_rec._name.substr(0, 4) == _source_pdb.substr(0, 4)) {
                    _records[l].first.clear(); _records[l].first.push_back(templ_rec);
                    num = 1;
                    break;
                }
            } else if (nuc2d::pure_ss(nuc2d::lower_ss(templ_rec._ss, 1)) == lower_ss) {
                templ_rec._score = (templ_rec._ss == ss ? 5 : 0);
                if (templ_rec._name.substr(0, 4) == _name.substr(0, 4)) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                if (templ_rec._family == family && family != "other") {
                    if (_is_test) continue; else templ_rec._score += 2;
                }
                for (int i = 0; i < templ_rec._seq.size(); i++) {
                    if (seq[i] == templ_rec._seq[i]) templ_rec._score++;
                }
                _records[l].first.push_back(templ_rec);
                num++;
                if (not _source_pdb.empty() and templ_rec._name.substr(0, 4) == _source_pdb.substr(0, 4)) {
                    _records[l].first.clear(); _records[l].first.push_back(templ_rec);
                    num = 1;
                    break;
                }
            }
        }
        ifile.close();
        if (_records[l].first.size() == 0) {
            _records[l].first.push_back();
        }
        std::sort(_records[l].first.begin(), _records[l].first.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
        _loop_nums_table.push_back(std::make_tuple(l, num_sons, num, 0));
    }

    void find_helix_records(nuc2d::loop *l) {
        if (l->s.empty()) return;

        std::string seq = l->s.seq(), ss = l->s.ss(), family = _family;
        int len = l->s.len();

        std::string info_file = _lib + "/" + _type + "/" + "records/helix";
        std::ifstream ifile(info_file.c_str());
        if (!ifile) throw "jian::nuc3d::FindTemplates error! Can't find '" + info_file + "'!";

        TemplRec templ_rec;
        int num = 0;
        while (ifile >> templ_rec._name >> templ_rec._len >> templ_rec._seq >> templ_rec._ss >> templ_rec._family) {
            if (templ_rec._len == len) {
                templ_rec._score = 0;
                if (templ_rec._name.substr(0, 4) == _name.substr(0, 4)) {
                    if (_is_test) continue; else templ_rec._score += 5;
                }
                if (templ_rec._family == family && family != "other") {
                    if (_is_test) continue; else templ_rec._score += 2;
                }
                for (int i = 0; i < templ_rec._seq.size(); i++) {
                    if (seq[i] == templ_rec._seq[i]) {
                        templ_rec._score++;
                    }
                }
                _records[l].second.push_back(std::move(templ_rec));
                num++;
            }
        }
        ifile.close();
        std::sort(_records[l].second.begin(), _records[l].second.end(), []( const TemplRec &loop1, const TemplRec &loop2) {
            return loop1._score > loop2._score; });
        for (auto &&tuple: _loop_nums_table) if (std::get<0>(tuple) == l) {std::get<3>(tuple) = num; break;}
    }

};

} // namespace nuc3d
} // namespace jian

#endif

