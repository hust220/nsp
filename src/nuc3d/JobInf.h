#ifndef JIAN_NUC3D_JOBINF_H
#define JIAN_NUC3D_JOBINF_H

#include <util/std.h>
#include <util/Par.h>
#include <nuc2d/N2D.h>
#include <nuc2d/util.h>

namespace jian {
namespace nuc3d {

class JobInf {
public:
    std::string _name = "assemble";
    std::string _seq;
    std::string _ss;
    std::string _lib;
    std::string _family = "other";
    std::string _type = "RNA";
    std::string _constraints;
    int _hinge = 2;
    int _num = 1;
    int _is_test = 0; // Is this a test case or not?
    std::string _method = "assemble";
    std::string _native;
    nuc2d::N2D _n2d;

    JobInf() {
        _lib = env("NSP");
    }

    JobInf(Par pars) {
        // ## Set sequence
        pars.count("sequence") && (_seq = pars["sequence"][0], 1);
        pars.count("seq") && (_seq = pars["seq"][0], 1);
        boost::to_upper(_seq);

        // ## Set secondary structure
        pars.count("secondary_structure") && (_ss = pars["secondary_structure"][0], 1);
        pars.count("ss") && (_ss = pars["ss"][0], 1);

        // ## Set library
        _lib = env("NSP"); // default library
        pars.count("library_path") && (_lib = pars["library_path"][0], 1);
        pars.count("lib") && (_lib = pars["lib"][0], 1);

        // ## Set job name
        pars.count("job_name") && (_name = pars["job_name"][0], 1);
        pars.count("job") && (_name = pars["job"][0], 1);

        // ## Set prediction number
        pars.count("number") && (_num = std::stoi(pars["number"][0]), 1);
        pars.count("num") && (_num = std::stoi(pars["num"][0]), 1);

        // ## Set hinge size
        if (pars.count("hinge")) _hinge = stoi(pars["hinge"][0]);

        // ## Set family
        pars.count("family") && (_family = pars["family"][0], 1);

        // ## Set prediction type
        pars.count("type") && (_type = pars["type"][0], 1);

        // ## Set constraints file
        pars.count("constraints") && (_constraints = pars["constraints"][0], 1);

        // ## Set meethod
        pars.count("method") && (_method = boost::to_lower_copy(pars["method"][0]), 1);

        // ## Set whether test
        pars.count("test") && (_is_test = ( boost::to_lower_copy(pars["test"][0]) == "yes" ? 1 : 0), 1);

        // ## Set the path of native pdb
        pars.count("native") && (_native = boost::to_lower_copy(pars["native"][0]), 1);

        // ## Check whether sequence and secondary are null
        _ss != "" || die("Please tell me the secondary structure!");
        _seq != "" || die("Please tell me the sequence!");

        // ## Check the length of the secondary structure and sequence
        if (!nuc2d::check_ss(_ss)) { 
            throw "The secondary structure includes illegal characters!";
        }
        if (nuc2d::len_ss(_ss) != _seq.size()) {
            throw "The length of the secondary structure and sequence should be equal!";
        }
    }

};

} // namespace nuc3d
} // namespace jian

#endif




