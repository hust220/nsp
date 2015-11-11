#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <nuc2d/util.h>
#include <pdb/util.h>
#include "Connect.h"

namespace jian {

namespace nuc3d {

class Assemble {
public:
    Assemble(Par);
    Assemble(string, string);
    void operator ()();

    string ss;
    string seq;
    int hinge_size = 2;
    string job = "assemble";
    string lib;
    string family = "other";
    string type = "RNA";
    Log log {"nsp.log"};
    nuc2d::N2D mol;
    string constraints;
    int num = 1;
    int is_test = 0; // Is this a test case or not?
    int view = 0;
    std::string _method = "assemble";

    /////////////////////////////////////
    /// Find templates
    int find_templates(nuc2d::loop *);
    vector<Model> find_loops(nuc2d::loop *, int = 1);
    vector<Model> find_helices(nuc2d::helix *, int = 1);
    void create_loops();

    int _max_loop_nums = 1000;
    map<nuc2d::loop *, pair<vector<Model>, vector<Model>>> templates;
    vector<tuple<nuc2d::loop *, int, int, int>> loop_nums_table;
    //////////////////////////////////////

    //////////////////////////////////////
    /// Assemble templates
    void ass_templates(int = 1);
    Model ass_templates(nuc2d::loop *);
    void rebuild_chains(Model &);
    //Model createHelix(nuc2d::helix *);
    //Model createHelix(string);

    map<nuc2d::loop *, pair<int, int>> templ;
    Connect connect;
    //////////////////////////////////////
};

} /// namespace nuc3d

} /// namespace jian

#endif




