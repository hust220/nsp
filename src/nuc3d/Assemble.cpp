#include "Assemble.h"
#include "BuildHelix.h"

namespace jian {

namespace nuc3d {

using namespace nuc2d;

Assemble::Assemble(Par pars) {
    /// set default parameter
    lib = env("NSP");

    /// read parameter file
    pars.count("sequence") && (seq = pars["sequence"][0], 1);
    pars.count("seq") && (seq = pars["seq"][0], 1);
    boost::to_upper(seq);

    pars.count("secondary_structure") && (ss = pars["secondary_structure"][0], 1);
    pars.count("ss") && (ss = pars["ss"][0], 1);

    pars.count("library_path") && (lib = pars["library_path"][0], 1);
    pars.count("lib") && (lib = pars["lib"][0], 1);

    pars.count("job_name") && (job = pars["job_name"][0], 1);
    pars.count("job") && (job = pars["job"][0], 1);

    pars.count("number") && (num = stoi(pars["number"][0]), 1);
    pars.count("num") && (num = stoi(pars["num"][0]), 1);

    if (pars.count("hinge")) {
        hinge_size = stoi(pars["hinge"][0]);
        connect._hinge_size = hinge_size;
    }

    pars.count("family") && (family = pars["family"][0], 1);
    pars.count("type") && (type = pars["type"][0], 1);
    pars.count("constraints") && (constraints = pars["constraints"][0], 1);
    pars.count("max_search_number") && (_max_loop_nums = stoi(pars["max_search_number"][0]), 1);
    pars.count("method") && (_method = boost::to_lower_copy(pars["method"][0]), 1);
    pars.count("test") && (is_test = (boost::to_lower_copy(pars["test"][0]) == "yes" ? 1 : 0), 1);
    pars.count("view") && (view = (boost::to_lower_copy(pars["view"][0]) == "yes" ? 1 : 0), 1);

    /// set library path
    if (upper(type) == "RNA") {
        lib += "/RNA/";
    } else if (upper(type) == "DNA") {
        lib += "/DNA/";
    } else if (lower(type) == "protein") {
        lib += "protein/";
    }

    ss != "" || die("Please tell me the secondary structure!");
    seq != "" || die("Please tell me the sequence!");

    /// check the length of the secondary structure and sequence
    count_if(begin(ss), end(ss), [](const char &c) {
        return c == '.' || c == '(' || c == ')' || c == '[' || c == ']';
    }) == seq.size() || die("The length of the secondary structure and sequence should be equal!");
}

Assemble::Assemble(string seq, string ss) {
    /// set default parameter
    lib = env("NSP");

    /// set secondary structure and sequence
    this->ss = ss;
    this->seq = seq;

    /// check the length of the secondary structure and sequence
    count_if(begin(ss), end(ss), [](const char &c) {
        return c == '.' || c == '(' || c == ')' || c == '[' || c == ']';
    }) == seq.size() || die("The length of the secondary structure and sequence should be equal!");
}

void Assemble::operator ()() {
    log("=========================================================");
    log("Job: " + job);
    log(Time::time() + "  start...");
    boost::timer t;
    log(seq);
    log(ss);
    log("----------------------------------------");

    log("Construct 2D structures...");
    mol.hinge_base_pair_num = hinge_size;
    mol(seq, ss);

    log("Find templates...");
    find_templates(mol.pseudo_head);

    log("Assemble templates...");
    ass_templates(num);

    log("----------------------------------------");
    log(Time::time() + "  done...");
    log("Time elapsed: " + std::to_string(t.elapsed()) + "s");
    log("=========================================================");
    log("\n\n");
}

} /// namespace nuc3d

} /// namespace jian


