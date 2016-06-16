#pragma once

#include "Residue.hpp"
#include <iostream>
#include <string>
#include <deque>

namespace jian {

class Chain : public std::deque<Residue> {
public:
    std::string name = "X";
    std::string type = "unknown";
    std::string model_name;

    Chain();

    template<typename T>
    Chain coarse_grained(T &&names) const {
        Chain c;
        c.name = name;
        for (auto &&res : *this) {
            c.push_back(res.coarse_grained(names));
        }
        return c;
    }

};

std::ostream &operator <<(std::ostream &output, const Chain &chain);
Chain residues_from_file(const std::string &file_name);
void residues_to_file(const Chain &chain, const std::string &file_name);
void append_chain_to_file(const Chain &chain, const std::string &file_name, int n);

} // namespace jian

