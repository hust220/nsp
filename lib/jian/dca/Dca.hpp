#pragma once

#include <utility>
#include <algorithm>
#include <list>
#include <tuple>
#include <map>
#include <deque>
#include <string>
#include <jian/utils/Par.hpp>
#include <jian/utils/file.hpp>
#include <jian/utils/log.hpp>
#include <jian/matrix.hpp>
#include "../utils/Factory.hpp"

namespace jian {
namespace dca {

class Dca {
public:
    using seqs_t = std::deque<std::string>;
    using creator_t = Dca *();

    int N, M, q;
	Mati align;
    Matf seqids, fi, Pi, DI;
    Vecf ma;
    Mat4 fij, Pij;
    float theta = 0.8f;
    float Meff;
    float pseudocount_weight = 0.5f;
    std::string m_seq;
    std::vector<int> m_indices;
    std::vector<char> m_symbols {'A', 'U', 'C', 'G'};
    std::map<char, int> m_map_symbols;

    ~Dca();
    void fastaread(std::string fastafile, seqs_t &seqs);
    void trim_seqs(seqs_t &seqs, int n);
    void init(std::string mol_type, std::string Rfamfile, int n);
    float seqid(int a, int b);
    void cal_seqids();
    void cal_ma();
    void cal_meff();
    void cal_fi();
    void cal_fij();
    void calculate_f();
    void calculate_P();
    void calculate_DI(std::string out_file);
    Dca &run(std::string mol_type, std::string input, std::string out_file, int n);
    virtual void calculate_eij() = 0;
    virtual float cal_di(int i, int j) = 0;
    virtual void set_step(float) {}
};

#define REG_DCA_FAC(name, type) REGISTER_FACTORY(Dca::creator_t, name, type)
using FacDca = Factory<Dca::creator_t>;

} // namespace dca
} // namespace jian

