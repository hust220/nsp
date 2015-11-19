#include "Assemble.h"
#include "Transform.h"
#include "AddPhos.h"

namespace jian {

namespace nuc3d {

using namespace nuc2d;

void Assemble::ass_templates(int n) {
    for (auto &temp_pair: templates) {
        templ[temp_pair.first] = {0, 0};
    }

    int flag = 0;
    int index = 0;
    while (flag == 0 && index < n) {
        log("Assemble model " + to_string(index + 1), 2);
        auto model = ass_templates(mol.pseudo_head);

        log("Rebuild chains...", 4);
        rebuild_chains(model);

        log("Mutate...", 4);
        model = Transform(model)(type, seq);

        log("Add phosphate group...", 4);
        AddPhos()(model);

        if (type == "DNA") {
            DNA(model).write(job + "-" + to_string(index + 1) + ".pdb");
        } else if (type == "RNA") {
            RNA(model).write(job + "-" + to_string(index + 1) + ".pdb");
        }

        flag = 1;
        for (auto &temp_pair: templates) {
            if (flag) {
                if (temp_pair.first->head != NULL) templ[temp_pair.first].first++;
                flag = 0;
                if (templ[temp_pair.first].first == templates[temp_pair.first].first.size()) {
                    templ[temp_pair.first].first = 0;
                    templ[temp_pair.first].second++;
                    if (templ[temp_pair.first].second == templates[temp_pair.first].second.size()) {
                        templ[temp_pair.first].second = 0;
                        flag = 1;
                    }
                }
            }
        }

        index++;
    }
}

Model Assemble::ass_templates(loop *l) {
    /// find templates of itself
    int a1 = l->s.getLen() - 1;
    int a2 = a1 + 1;
    Model model, helix_model, loop_model;
    if (l->s.head == NULL) {
        model = templates[l].first[templ[l].first];
    } else {
        helix_model = templates[l].second[templ[l].second];
        if (l->head == NULL) {
            return helix_model;
        } else {
            loop_model = templates[l].first[templ[l].first];
            model = connect(helix_model, loop_model, a1, a2);
        }
    }

    /// assemble
    int n = (l->s.head == NULL) ? 0 : (l->s.getLen() - mol.hinge_base_pair_num);
    int hinge_index = 0;
    for (loop *tempLoop = l->son; tempLoop != NULL; tempLoop = tempLoop->brother, hinge_index++) {
        auto temp_model = ass_templates(tempLoop);
        model = connect(model, temp_model, l->hinges[hinge_index].first + n, l->hinges[hinge_index].second + n);
        n += temp_model.res_nums() - 4;
    }
    return model;
}

//Model Assemble::createHelix(helix *s) {
//    string seq1 = "";
//    string seq2 = "";
//    for (bp *b = s->head; b != NULL; b = b->next) {
//        seq1 += b->res1.name;
//        seq2 += b->res2.name;
//    }
//    for (int i = seq2.size() - 1; i >= 0; i--) {
//        seq1 += seq2[i];
//    }
//    return createHelix(seq1);
//}
//
//Model Assemble::createHelix(string seq) {
//    std::stringstream strstr;
//    std::string file_name;
//
//    if (seq.size() < 4) {
//        std::cerr << "Assemble::createHelix error! The length of the helix"
//                  << " to be create (" << seq << ")should not be less than 4!" << std::endl;
//        exit(1);
//    } else if (seq.size() % 2 == 1) {
//        std::cerr << "Assemble::createHelix error! The length of the helix"
//                  << " to be create (" << seq << ")should be even!" << std::endl;
//        exit(1);
//    } else if (seq.size() == 4 || seq.size() == 6) {
//        strstr << lib << "/basepair/" << seq << ".pdb";
//        strstr >> file_name;
//        ifstream ifile(file_name.c_str());
//        if (!ifile) {
//            //strstr.str("");
//            strstr.clear();
//            strstr.seekp(0);
//            strstr.seekg(0);
//            strstr << lib << "/basepair/XXXXXX.pdb";
//            strstr >> file_name;
//        }
//        ifile.close();
//        return new Model(file_name);
//    } else {
//        strstr << lib << "/basepair/" << seq.substr(0, 3) 
//               << seq.substr(seq.size() - 3, 3) << ".pdb";
//        strstr >> file_name;
//        ifstream ifile(file_name.c_str());
//        if (!ifile) {
//            //strstr.str("");
//            strstr.clear();
//            strstr.seekp(0);
//            strstr.seekg(0);
//            strstr << lib << "/basepair/XXXXXX.pdb";
//            strstr >> file_name;
//        }
//        ifile.close();
//        return connect(Model(file_name), createHelix(seq.substr(1, seq.size() - 2)), 2, 3);
//    }
//}
//
void Assemble::rebuild_chains(Model &model) {
    /// Construct new chains
    vector<Chain> temp_chains;
    Chain temp_chain;
    int res_num = 0;
    for (auto &chain: model.chains) {
        for (auto &residue: chain.residues) {
            if (ss[res_num] == '&') {
                res_num++;
                temp_chains.push_back(temp_chain);
                temp_chain.residues.clear();
            } else if (ss[res_num] == ')' && ss[res_num - 1] == '(' && res_num > 0) {
                temp_chains.push_back(temp_chain);
                temp_chain.residues.clear();
            }
            temp_chain.residues.push_back(residue);
            res_num++;
        }
    }
    temp_chains.push_back(temp_chain);
    model.chains = temp_chains;
}

} /// namespace nuc3d

} /// namespace jian



