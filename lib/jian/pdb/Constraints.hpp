#pragma once

#include "../utils/ls.hpp"

namespace jian {

struct Constraint {
    std::vector<int> key;
    double value;
    double min, max;
    double weight;
};

Constraint make_contact(int i, int j, double weight = 1);
Constraint make_distance(int i, int j, double value, double weight = 1);
Constraint make_angle(int i, int j, int k, double value, double weight = 1);
Constraint make_dihedral(int i, int j, int k, int l, double value, double weight = 1);

class Constraints {
public:
    std::deque<Constraint> contacts;
    std::deque<Constraint> distances;
    std::deque<Constraint> angles;
    std::deque<Constraint> dihedrals;

    void add_contact(const Constraint &ct);
    void add_distance(const Constraint &ct);
    void add_angle(const Constraint &ct);
    void add_dihedral(const Constraint &ct);
    bool has_contact(int i, int j);
    bool has_distance(int i, int j);
    bool has_angle(int i, int j, int k);
    bool has_dihedral(int i, int j, int k, int l);
    void read_distances_file(const std::string &f);

};

Constraints read_constraints(const std::string &f);

} // namespace jian

