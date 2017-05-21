#include "rtsp_templ_rec.hpp"
#include "rna_mc_templates.hpp"

BEGIN_JN

static void find_loop_records(SSE *l, Deque<TemplRec> &records, Int max_num = -1, const Deque<Str> &sources = {}, const Deque<Str> &excludes = {})
{
    if (l->loop.empty()) return;

    Str seq = l->loop.seq(), ss = l->loop.ss(), p_ss = NASS::pure_ss(ss), lower_ss = NASS::lower_ss(p_ss, 1);

    int num_sons = l->num_sons();

    STD_ ostringstream stream;
    stream << Env::lib() << "/RNA/records/" << (l->is_open() ? "open_" : "") << "loop";
    STD_ ifstream ifile;
    FOPEN(ifile, stream.str());

    TemplRec templ_rec;
    S line;
    int num = 0;
    while (std::getline(ifile, line) && set_loop_rec(templ_rec, line)) {
        if (std::find(excludes.begin(), excludes.end(), templ_rec.src) != excludes.end()) {
            continue;
        } else if ((!sources.empty()) && std::find(sources.begin(), sources.end(), templ_rec.src) == sources.end()) {
            continue;
        } else if (NASS::pure_ss(NASS::lower_ss(templ_rec.ss, 1)) == lower_ss) {
            templ_rec.score = (templ_rec.ss == ss ? 5 : 0);
            for (uint i = 0; i < templ_rec.seq.size(); i++) {
                if (seq[i] == templ_rec.seq[i]) {
                    if (ss[i] == templ_rec.ss[i] && ss[i] != '.' && ss[i] != '(' && ss[i] != ')') {
                        templ_rec.score += 2;
                    } else if (ss[i] == '(' || ss[i] == ')') {
                        templ_rec.score += 0.2;
                    } else {
                        templ_rec.score += 1;
                    }
                }
            }
            records.push_back(templ_rec);
            num++;
        }
    }
    ifile.close();
    std::sort(records.begin(), records.end(), [](auto &&r1, auto &&r2) { return r1.score > r2.score; });
    if (max_num != -1 && size(records) > max_num) records.erase(std::next(records.begin(), max_num+1), records.end());
}

static void find_helix_records(SSE *l, Deque<TemplRec> &records, Int max_num) {
    if (l->helix.empty()) return;

    S seq = l->helix.seq();
    S ss = l->helix.ss();

    int len = size(l->helix);

    STD_ ostringstream stream;
    stream << Env::lib() << "/RNA/records/helix";
    STD_ ifstream ifile;
    FOPEN(ifile, stream.str());

    TemplRec templ_rec;
    Str line;
    int num = 0;
    while (std::getline(ifile, line) && set_helix_rec(templ_rec, line)) {
        if (templ_rec.len == len) {
            templ_rec.score = 0;
            for (uint i = 0; i < templ_rec.seq.size(); i++) {
                if (seq[i] == templ_rec.seq[i]) {
                    templ_rec.score++;
                }
            }
            records.push_back(std::move(templ_rec));
            num++;
        }
    }
    ifile.close();

    std::sort(records.begin(), records.end(), [](auto &&r1, auto &&r2) { return r1.score > r2.score; });
    if (max_num != -1 && size(records) > max_num) records.erase(std::next(records.begin(), max_num+1), records.end());
}

void dhmc_set_templates(DHMC &m) {
    for_ (sse, tree_nodes(m.m_sst.head)) {
//        std::cout << *sse << std::endl;
//        std::cout << "set loop records" << std::endl;
        find_loop_records(sse, m.m_loop_records[sse], 50, m.m_sources, m.m_excludes);
//        std::cout << "set helix records" << std::endl;
        find_helix_records(sse, m.m_helix_records[sse], 50);
    } end_;
}

void dhmc_print_templates(DHMC &m) {
    m.log << "# Print templates..." << std::endl;
    for_ (sse, tree_nodes(m.m_sst.head)) {
        m.log << *sse << std::endl;
        m.log << "helix records: " << join(" ", map_s<Deque<Str>>([](auto && r){return r.name;}, m.m_helix_records[sse])) << std::endl;
        m.log << "loop records: " << join(" ", map_s<Deque<Str>>([](auto && r){return r.name;}, m.m_loop_records[sse])) << std::endl;
    } end_;
}

END_JN

