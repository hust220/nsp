#include <jian/utils/string.hpp>
#include <jian/utils/Env.hpp>
#include <jian/utils/file.hpp>
#include "TemplRec.hpp"

BEGIN_JN

static void rec_set_src(TemplRec &rec) {
    int pos = rec.name.find_first_of("-");
    if (pos == Str::npos) {
        rec.src = jian::upper(rec.name);
    } else {
        rec.src = jian::upper(rec.name.substr(0, pos));
    }
}

Bool set_loop_rec(TemplRec &rec, const Str &s) {
    tokenize_v v;
    jian::tokenize(s, v, " ");
    if (v.size() != 5) {
        return false;
    } else {
        rec.name = v[0];
        rec_set_src(rec);
        rec.type = JN_INT(v[1]);
        rec.seq = v[2];
        rec.ss = v[3];
        rec.family = v[4];
        return true;
    }
}

Bool set_helix_rec(TemplRec &rec, const Str &s) {
    tokenize_v v;
    jian::tokenize(s, v, " ");
    if (v.size() != 5) {
        return false;
    } else {
        rec.name = v[0];
        rec_set_src(rec);
        rec.len = JN_INT(v[1]);
        rec.seq = v[2];
        rec.ss = v[3];
        rec.family = v[4];
        return true;
    }
}

static Str records_path() {
    return to_str(Env::lib(), "/RNA/records/");
}

Deque<TemplRec> find_loop_templates(Str seq, Str ss, Str src) {
    if (seq.empty()) return {};

    Str p_ss = NASS::pure_ss(ss), lower_ss = NASS::lower_ss(p_ss, 1);

    Str info_file = to_str(records_path(), "loop");
    std::ifstream ifile;
    FOPEN(ifile, info_file);

    Deque<TemplRec> records;
    TemplRec rec;
    Str line;
    while (std::getline(ifile, line) && set_loop_rec(rec, line)) {
        if (NASS::pure_ss(NASS::lower_ss(rec.ss, 1)) == lower_ss) {
            rec.score = (rec.ss == ss ? 5 : 0);
            rec.score += (rec.src == src ? 5 : 0);
            for (Int i = 0; i < size(rec.seq); i++) {
                if (seq[i] == rec.seq[i]) {
                    if (ss[i] == rec.ss[i] && ss[i] != '.' && ss[i] != '(' && ss[i] != ')') {
                        rec.score += 2;
                    } else if (ss[i] == '(' || ss[i] == ')') {
                        rec.score += 0.2;
                    } else {
                        rec.score += 1;
                    }
                }
            }
            records.push_back(rec);
        }
    }
    ifile.close();
    std::sort(records.begin(), records.end(), [](auto && r1, auto && r2) {
            return r1.score > r2.score; });
    return records;
}

static Str helix_ss(Int n) {
    std::stringstream stream;
    for (Int i = 0; i < n/2; i++) stream << '(';
    for (Int i = 0; i < n/2; i++) stream << ')';
    return stream.str();
}

Deque<TemplRec> find_helix_templates(Str seq, Str src) {
    Int l = size(seq);

    if (l < 4) return {};

    Str ss = helix_ss(l);

    Str info_file = to_str(records_path(), "helix");;
    std::ifstream ifile(info_file.c_str());
    if (!ifile) throw "jian::FindTemplates error! Can't find '" + info_file + "'!";

    Deque<TemplRec> records;
    TemplRec rec;
    Str line;
    while (std::getline(ifile, line) && set_helix_rec(rec, line)) {
        if (rec.len == l) {
            rec.score = (rec.src == src ? 5 : 0);
            for (Int i = 0; i < size(rec.seq); i++) {
                if (seq[i] == rec.seq[i]) {
                    rec.score++;
                }
            }
            records.push_back(rec);
        }
    }
    ifile.close();
    std::sort(records.begin(), records.end(), [](auto && r1, auto && r2) {
            return r1.score > r2.score; });
    return records;
}

END_JN

