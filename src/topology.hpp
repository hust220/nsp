#pragma once

namespace jian {

class Topology {
public:
    double min_inter_distance;

    /**
     * Get the number of atoms given the sequence.
     */
    virtual int get_num_atoms(const std::string &seq) const = 0;
};

class Topology1p : Topology {
public:
    Topology1p() : Topology() {
        this->min_inter_distance = 6.1;
    }

    /**
     * Get the number of atoms given the sequence.
     */
    virtual int get_num_atoms(const std::string &seq) const {
        return seq.size();
    }
};

Topology *make_topology(const std::string &type) {
    if (type == "1p") {
        return new Topology();
    }
    else {
        throw "error";
    }
}

}

