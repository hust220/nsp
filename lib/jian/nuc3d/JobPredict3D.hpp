#pragma once

#include <string>
#include <ctime>
#include "../pdb/Constraints.hpp"
#include "../nuc2d/SSTree.hpp"

namespace jian {

class Par;

class JobPredict3D {
public:
    std::string _name = "assemble";
    std::string _seq;
    std::string _ss;
    std::string _lib;
    std::string _family = "other";
    std::string _type = "RNA";
    std::string _file_constraints;
    Constraints _constraints;
    std::deque<std::string> m_disused_pdbs;
    std::time_t _start_time, _end_time;
    int _hinge = 2;
    int _num = 1;
    int _num_sampling = 100;
    bool _is_test = false; // Is this a test case or not?
    std::string _method = "FA";
    std::string _native;
    SSTree _ss_tree;
    std::string _strategy = "loose";
    std::string _source_pdb;
    double m_seed = 11;
//    std::string m_out;
//    bool m_no_mc = false;
    bool m_sample_hairpin = false;
    bool m_sample_il = false;
    std::string m_mode = "auto";
    const Par *_par = NULL;
    std::string m_cmd;

    JobPredict3D();
    JobPredict3D(const Par &pars);
    void set_constraints();
    void display_start_information();
    void display_end_information();

};

} // namespace jian

