#include "Assemble.h"
#include "LoopModelling2.h"

namespace jian {

namespace nuc3d {

using namespace nuc2d;

int Assemble::find_templates(loop *l) {
    if (l == NULL) {
        return 0;
    }

    /// find templates of helix
    if (l->head != NULL) {
        templates[l].first = find_loops(l, _max_loop_nums);
    }
    
    /// find templates of loop
    if (l->s.head != NULL) {
        if (l == mol.pseudo_head && l->head == NULL) {
            templates[l].second = find_helices(&(l->s), num);
        } else {
            templates[l].second = find_helices(&(l->s), 1);
        }
    }
    loop_nums_table.push_back(make_tuple(l, l->hinges.size(), templates[l].first.size(), templates[l].second.size()));

    /// find templates of son
    find_templates(l->son);

    /// find templates of brother
    find_templates(l->brother);

    if (l == mol.pseudo_head) {
        sort(begin(loop_nums_table), end(loop_nums_table), [](
            const tuple<loop *, int, int, int> &tuple1, const tuple<loop *, int, int, int> &tuple2
        ){
            return get<1>(tuple1) > get<1>(tuple2);
        });
        auto product = [&] {
            return accumulate(begin(loop_nums_table), end(loop_nums_table), 1, [&](
                int a, const tuple<loop *, int, int, int> &t
            ){
                if (get<0>(t)->head != NULL) {
                    return (a > num ? num: a * get<2>(t));
                } else {
                    return a;
                }
            });
        };
        for (auto &t: loop_nums_table) {
            if (get<0>(t)->head != NULL && get<2>(t) == 0) get<2>(t)++;
            if (get<0>(t)->s.head != NULL && get<3>(t) == 0) get<3>(t)++;
        }
        if (l->head != NULL && product() < num) {
            while (true) {
                for (auto &t: loop_nums_table) {
                    if (get<0>(t)->head != NULL) {
                        get<2>(t)++;
                    } else {
                        continue;
                    }
                    if (product() >= num) goto endwhile;
                }
            } endwhile:;
        }
        create_loops();
    }
}

vector<Model> Assemble::find_loops(loop *l, int nums) {
    /// select method
    if (_method == "dg")
        return std::vector<Model>();

    /// find loop from loop library
    string l_seq = l->getSeq();
    string l_ss = l->getSS();
    string l_purified_ss;
    copy_if(l_ss.begin(), l_ss.end(), back_inserter(l_purified_ss), [](char c) {
        return c != '&';
    });
    string l_family = family;
    int hinge_num = l->hinges.size();

    string info_file = lib;
    if (count_if(l_ss.begin(), l_ss.end(), [](char c) {
        return c == '[' || c == ']';
    })) {
        info_file += "/info-" + std::to_string(hinge_size) + "/pseudo_knot";
    } else {
        info_file += "/info-" + std::to_string(hinge_size) + "/loop-" + to_string(hinge_num);
    }

    ifstream ifile(info_file.c_str());
    ifile || die(string("Assemble::findLoop error! Can't find '") + info_file + "'!");

    class LoopItem {
    public:
        string name;
        string seq;
        string ss;
        string family;
        int score;
    };
    LoopItem loop_item;

    vector<LoopItem> loop_sets;
    while (ifile >> loop_item.name >> loop_item.seq >> loop_item.ss >> loop_item.family) {
        if (loop_item.ss == l_ss) {
            loop_item.score = 0;

            // If this is a test case, then skip the templates from itself;
            // If it's not a test case, then increase the score by 2.
            if (loop_item.name.substr(5).compare(job.substr(0, 4)) == 0) {
                if (is_test) {
                    continue;
                } else {
                    loop_item.score += 5;
                }
            }

            // If this is a test case, then skip the templates from 
            // Models the family of which is the same;
            // If it's not a test case, then increase the score by 2.
            if (l_family == loop_item.family && l_family != "other") {
                if (is_test) {
                    continue;
                } else {
                    loop_item.score += 2;
                }
            }

            for (int i = 0; i < l_purified_ss.size(); i++) {
                if (l_seq[i] == loop_item.seq[i]) {
                    if (set<char>{'(', ')', '[', ']'}.count(l_purified_ss[i])) {
                        loop_item.score += 0.5;
                    } else {
                        loop_item.score += 1;
                    }
                }
            }

            loop_sets.push_back(loop_item);
            if (loop_sets.size() >= nums) break;
        }
    }
    ifile.close();
    sort(loop_sets.begin(), loop_sets.end(), [](
        const LoopItem &loop1, const LoopItem &loop2) {
            return loop1.score > loop2.score;
        }
    );

    log("for loop: " + l_seq + " " + l_ss, 2);

    vector<Model> loop_models;
    for (int i = 0; i < loop_sets.size(); i++) {
        // get loop name
        string file_name = lib + "/loop-" + std::to_string(hinge_size) + "/" + loop_sets[i].name + ".pdb";

        // log
        log("find: " + loop_sets[i].name, 4);
        log("sequence: " + loop_sets[i].seq, 6);
        log("secondary structure: " + loop_sets[i].ss, 6);
        log("family: " + loop_sets[i].family, 6);
        log("score: " + to_string(loop_sets[i].score), 6);

        loop_models.push_back(Model(file_name));
    }
    if (loop_sets.empty()) {
        log("find no appropriate templates!", 4);
    }

    return loop_models;
}

vector<Model> Assemble::find_helices(helix *s, int nums) {
    bp *b;
    int i;
    string seq1, seq2;

    // find helix info
    int len = s->getLen();
    string name = lib;
    name += "info-" + std::to_string(hinge_size) + "/helix";
    ifstream ifile(name.c_str());
    helixInfo *hi = new helixInfo;
    helixInfoList *hiList = new helixInfoList;

    class HelixItem {
    public:
        string name;
        int len;
        string seq;
        string ss;
        string src;
        int score;
    };
    HelixItem helix_item;

    vector<HelixItem> helix_sets;
    while (ifile >> helix_item.name >> helix_item.len >> helix_item.seq >> helix_item.ss >> helix_item.src) {
    //while (ifile >> hi->name >> hi->n >> hi->seq >> hi->ss >> hi->src) {
        if (len == helix_item.len) {
            int temp;
            double score;
            for (b = s->head, temp = 0, score = 0; b != NULL; b = b->next) {
                if (b->res1.name == helix_item.seq[temp]) {
                    if (temp == 0 || temp == helix_item.len - 1) {
                        score += 1 + 1 / (double) (2 * helix_item.len);
                    } else {
                        score += 1 / (double) (2 * helix_item.len);
                    }
                }
                if (b->res2.name == helix_item.seq[2 * len - 1 - temp]) {
                    if (temp == 0 || temp == helix_item.len - 1) {
                        score += 1 + 1 / (double) (2 * helix_item.len);
                    } else {
                        score += 1 / (double) (2 * helix_item.len);
                    }
                }
                temp++;
            }
            helix_item.score = score;
            helix_sets.push_back(helix_item);
            if (helix_sets.size() == nums) break;

//            hiList->add(hi, score);
//            hi = new helixInfo;
//            if (hiList->getLen() == nums) {
//                break;
//            }
        }
//        delete hi;
//        hi = new helixInfo;
    }
    ifile.close();

    log("for helix: " + s->seq() + " " + s->ss(), 2);

//    if (helix_sets.empty()) {
//    //if (hiList->head == NULL) {
//        log("create helix...", 4);
//        seq1 = "";
//        seq2 = "";
//        for (b = s->head; b != NULL; b = b->next) {
//            seq1 += b->res1.name;
//            seq2 += b->res2.name;
//        }
//        for (i = seq2.size() - 1; i >= 0; i--) {
//            seq1 += seq2[i];
//        }
//        return vector<Model>{createHelix(seq1)};
//    } else {
//        hi = hiList->head->hi;
//        name = lib + "helix/" + hi->name + ".pdb";
//        log("find: " + name, 4);
//        log("sequence: " + hi->seq, 6);
//        log("score: " + to_string(hiList->head->score), 6);
//        return vector<Model>{Model(name)};
//    }

    vector<Model> helices;
    for (int i = 0; i < helix_sets.size(); i++) {
        // get loop name
        string file_name = lib + "/helix-" + std::to_string(hinge_size) + "/" + helix_sets[i].name + ".pdb";

        // log
        log("find: " + helix_sets[i].name, 4);
        log("sequence: " + helix_sets[i].seq, 6);
        log("source: " + helix_sets[i].src, 6);
        log("score: " + to_string(helix_sets[i].score), 6);

        helices.push_back(Model(file_name));
    }
    if (helix_sets.empty()) {
        log("find no appropriate templates!", 4);
    }

    return helices;
}

void Assemble::create_loops() {
    log("Create templates...");
    for (auto &t: loop_nums_table) {
        loop *l; int l_type; int l_nums; int h_nums;
        tie(l, l_type, l_nums, h_nums) = t;
        if (l->head != NULL) {
            int diff = l_nums - templates[l].first.size();
            if (diff > 0) {
                string l_seq = l->getSeq();
                string l_ss = l->getSS();
                log("for loop: " + l_seq + " " + l_ss, 2);
                LoopModelling2 lm(type);
                lm._view = view;
                for (int i = 0; i < diff; i++) {
                    templates[l].first.push_back(Model(lm(l_seq, l_ss, constraints)));
                    log("create a loop template (penalty energy: " + to_string(lm.dg.E) + ").", 4);
                }
            }
        }
        if (l->s.head != NULL) {
            int diff = h_nums - templates[l].second.size();
            if (diff > 0) {
                string h_seq = l->s.seq();
                string h_ss = l->s.ss();
                log("for helix: " + h_seq + " " + h_ss, 2);
                LoopModelling2 lm(type);
                lm._view = view;
                for (int i = 0; i < diff; i++) {
                    templates[l].second.push_back(Model(lm(h_seq, h_ss)));
                    log("create a helix template (penalty energy: " + to_string(lm.dg.E) + ").", 4);
                }
            }
        }
    }
}

} /// namespace nuc3d

} /// namespace jian

