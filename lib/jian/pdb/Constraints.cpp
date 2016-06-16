#include "Constraints.hpp"
#include "../utils/file.hpp"
#include "../utils/ls.hpp"

namespace jian {

template<typename T>
static bool equal(T &&ls, int i,  int j) {
    return ls.size() == 2 && (ls[0] == i && ls[1] == j || ls[0] == j && ls[1] == i);
}

template<typename T>
static bool equal(T &&ls, int i,  int j, int k) {
    return ls.size() == 3 && ls[1] == j && (ls[0] == i && ls[2] == k || ls[0] == k && ls[2] == i);
}

template<typename T>
static bool equal(T &&ls, int i,  int j, int k,  int l) {
    return ls.size() == 4 && (ls[0] == i && ls[3] == l || ls[0] == l && ls[3] == i) && 
                             (ls[1] == j && ls[2] == k || ls[1] == k && ls[2] == j);
}

Constraint make_contact(int i, int j, double weight) {
    Constraint c; append(c.key, i, j); c.weight = weight; return c;
}

Constraint make_distance(int i, int j, double value, double weight) {
    Constraint c; append(c.key, i, j); c.value = value; c.weight = weight; return c;
}

Constraint make_angle(int i, int j, int k, double value, double weight) {
    Constraint c; append(c.key, i, j, k); c.value = value; c.weight = weight; return c;
}

Constraint make_dihedral(int i, int j, int k, int l, double value, double weight) {
    Constraint c; append(c.key, i, j, k, l); c.value = value; c.weight = weight; return c;
}

void Constraints::add_contact(const Constraint &ct) {
    contacts.push_back(ct);
}

void Constraints::add_distance(const Constraint &ct) {
    distances.push_back(ct);
}

void Constraints::add_angle(const Constraint &ct) {
    angles.push_back(ct);
}

void Constraints::add_dihedral(const Constraint &ct) {
    dihedrals.push_back(ct);
}

bool Constraints::has_contact(int i, int j) {
    for (auto &&c : contacts) if (equal(c.key, i, j)) return true; return false;
}

bool Constraints::has_distance(int i, int j) {
    for (auto &&c : distances) if (equal(c.key, i, j)) return true; return false;
}

bool Constraints::has_angle(int i, int j, int k) {
    for (auto &&c : angles) if (equal(c.key, i, j, k)) return true; return false;
}

bool Constraints::has_dihedral(int i, int j, int k, int l) {
    for (auto &&c : dihedrals) if (equal(c.key, i, j, k, l)) return true; return false;
}

void Constraints::read_distances_file(const std::string &f) {
    if (! f.empty()) {
        EACH_SPLIT_LINE(f.c_str(), " ", 
            if (F.size() == 3) add_distance(make_distance(JN_INT(F[0])-1, JN_INT(F[1])-1, JN_DBL(F[2])));
        );
    }
}

Constraints read_constraints(const std::string &f) {
    Constraints c; if (f == "") return c;
    EACH_SPLIT_LINE(f.c_str(), " ", 
        if (F.size() == 3) c.add_distance(make_distance(JN_INT(F[0])-1, JN_INT(F[1])-1, JN_DBL(F[2])));
    );
    return c;
}

} // namespace jian

