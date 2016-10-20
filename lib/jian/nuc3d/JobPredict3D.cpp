#include <sstream>
#include "JobPredict3D.hpp"
#include "../utils/Debug.hpp"
#include "../utils/Par.hpp"
#include "../utils/Env.hpp"
#include "../utils/rand.hpp"
#include "../utils/log.hpp"
#include "../utils/string.hpp"

namespace jian {

void JobPredict3D::init(const Par &pars) {
    _par = new Par(pars);
    _lib = Env::lib();
    m_cmd = (*_par)[1];
    pars.set(_seq, "seq", "sequence");
    jian::to_upper(_seq);
    LOG << "# Reading secondary structure" << std::endl;
    read_ss();
    pars.set(_lib, "lib", "library_path");
    pars.set(_name, "job_name", "job", "name");
    pars.set(_num, "n", "num", "number");
    pars.set(_hinge, "h", "hinge");
    pars.set(_strategy, "strategy");
    pars.set(_family, "family");
    pars.set(_type, "t", "type", "mol_type");
    pars.set(_file_constraints, "c", "constraints");
    pars.set(m_seed, "seed");
    pars.set(m_mode, "mode");
    seed(m_seed);
//    pars.set(m_out, "out");
    set_constraints();
    if (pars.has("disused_pdbs")) {
        m_disused_pdbs = pars["disused_pdbs"];
    }
    for (auto && s : m_disused_pdbs) jian::to_upper(s);
    pars.set(_method, "method");
    pars.set(_is_test, "test");
    pars.set(_native, "native");
    pars.set(_source_pdb, "source_pdb");
    jian::to_upper(_source_pdb);
    if (pars.has("sample_hp")) m_sample_hairpin = true;
    if (pars.has("sample_il")) m_sample_il = true;

    if (_ss == "") throw "Please tell me the secondary structure!";
    if (_seq == "") throw "Please tell me the sequence!";

}

JobPredict3D::~JobPredict3D() {
    delete _par;
}

void JobPredict3D::set_constraints() {
    if (! _file_constraints.empty()) {
        _constraints.read_distances_file(_file_constraints);
    }
}

void JobPredict3D::display_start_information() {
    _start_time = std::time(nullptr);
    LOG << "=========================================================\n"
              << "New Job: " << _name << '\n'
              << "Time: " << std::asctime(std::localtime(&(_start_time))) << '\n'
              << "Sequence: " << _seq << '\n'
              << "2D Structure: " << _ss << '\n'
              << "Molecular Type: " << _type << '\n'
              << "Number: " << _num << '\n'
              << "Method: " << _method << '\n'
              << "----------------------------------------\n";
}

void JobPredict3D::display_end_information() {
    _end_time = std::time(nullptr);
    LOG << "----------------------------------------\n"
              << "Finish Time: " << std::asctime(std::localtime(&(_end_time))) << '\n'
              << "Time elapsed: " << _end_time - _start_time << "s\n"
              << "=========================================================\n\n";
}

void JobPredict3D::read_ss() {
    _par->set(_ss, "ss", "secondary_structure");
}

} // namespace jian

