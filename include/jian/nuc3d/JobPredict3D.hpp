#pragma once

#include "../etl.hpp"
#include "../nuc2d.hpp"
#include "Env.hpp"

namespace jian {

class JobPredict3D {
public:
    std::string _name = "assemble";
    std::string _seq;
    std::string _ss;
    std::string _lib;
    std::string _family = "other";
    std::string _type = "RNA";
    std::string _constraints;
    std::string _disused_pdb;
    std::time_t _start_time, _end_time;
    int _hinge = 2;
    int _num = 1;
    int _num_sampling = 1000;
    bool _is_test = false; // Is this a test case or not?
    std::string _method = "FA";
    std::string _native;
    SSTree _ss_tree;
    std::string _strategy = "loose";
    std::string _source_pdb;

    JobPredict3D() : _lib(Env::lib()) {}

    JobPredict3D(Par pars) : _lib(Env::lib()) {
        pars.set(_seq, "seq", "sequence"); boost::to_upper(_seq);
        pars.set(_ss, "ss", "secondary_structure");
        pars.set(_lib, "lib", "library_path");
        pars.set(_name, "job_name", "job", "name");
        pars.set(_num, "n", "num", "number");
        pars.set(_hinge, "h", "hinge");
        pars.set(_strategy, "strategy");
        pars.set(_family, "family");
        pars.set(_type, "t", "type");
        pars.set(_constraints, "c", "constraints");
        pars.set(_disused_pdb, "disused_pdb");
        pars.set(_method, "method");
        pars.set(_is_test, "test");
        pars.set(_native, "native");
        pars.set(_source_pdb, "source_pdb");

        _ss != "" or die("Please tell me the secondary structure!");
        _seq != "" or die("Please tell me the sequence!");

        if (NucSS::len_ss(_ss) != _seq.size()) throw "The length of the secondary structure and sequence should be equal!";
    }

};

template<typename T>
void display_start_information_job(T &&job) {
    job._start_time = std::time(nullptr);
    Trace::log("=========================================================\n",
               "New Job: ", job._name, '\n',
               "Time: ", std::asctime(std::localtime(&(job._start_time))), '\n',
               "Sequence: ", job._seq, '\n',
               "2D Structure: ", job._ss, '\n',
               "Molecular Type: ", job._type, '\n',
               "Number: ", job._num, '\n',
               "Number of Sampling: ", job._num_sampling, '\n',
               "Method: ", job._method, '\n',
               "----------------------------------------\n");
}

template<typename T>
void display_end_information_job(T &&job) {
    job._end_time = std::time(nullptr);
    Trace::log("----------------------------------------\n",
               "Finish Time: ", std::asctime(std::localtime(&(job._end_time))), '\n',
               "Time elapsed: ", job._end_time - job._start_time, "s\n",
               "=========================================================\n\n");
}

} // namespace jian

