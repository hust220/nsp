#ifndef JIAN_EXTRACT_FRAGMENT_H
#define JIAN_EXTRACT_FRAGMENT_H

#include <jian/etl.h>
#include <jian/geom.h>
#include <jian/pdb.h>

namespace jian {

void extract_fragment(const jian::Model &model, int len) {
    int num_residues = jian::num_residues(model);
    int num_residue = 0;
    std::deque<jian::Point> deque1;
    std::deque<jian::Point> deque2;
    std::deque<jian::Point> deque3;
    std::deque<jian::Residue> deque;
    int file_index = 0;
    for (auto &&chain: model) {
        for (auto &&residue: chain) {
            deque.push_back(residue);
            deque1.push_back(pos(atom(residue, "O5*")));
            deque2.push_back(pos(atom(residue, "C4*")));
            deque3.push_back(pos(atom(residue, "O3*")));
            if (deque2.size() >= 2) {
                double dist = jian::geom::distance(deque1.back(), deque3[deque3.size() - 2]);
                if (dist > 4) {
                    deque1.erase(deque1.begin(), std::next(deque1.begin(), deque1.size() - 1));
                    deque2.erase(deque2.begin(), std::next(deque2.begin(), deque2.size() - 1));
                    deque3.erase(deque3.begin(), std::next(deque3.begin(), deque3.size() - 1));
                    deque.erase(deque.begin(), std::next(deque.begin(), deque.size() - 1));
                }
            }
            if (deque2.size() >= len) {
                file_index++;
                std::string file_name = model.name + "-frag-" + boost::lexical_cast<std::string>(len) + "-" + boost::lexical_cast<std::string>(file_index);
                std::cout << file_name << " ";
                for (int i = 0; i < len; i++) {
                    for (int j = i + 1; j < len; j++) {
                        std::cout << jian::geom::distance(deque2[i], deque2[j]) << ' ';
                    }
                }
                jian::Model new_model;
                new_model.push_back(jian::Chain());
                std::copy(deque.begin(), deque.end(), std::back_inserter(new_model[0]));
                write_pdb(new_model, file_name + ".pdb");
                std::cout << std::endl;
                deque1.pop_front();
                deque2.pop_front();
                deque3.pop_front();
                deque.pop_front();
            }
            num_residue++;
        }
    }
}

} // namespace jian

#endif
