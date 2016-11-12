#pragma once

#include "Dca.hpp"

namespace jian {
namespace dca {

class MpDca : public Dca {
public:
    Matf hi;
    Matf pi;
    Mat3 mij;
    Mat4 eij;
    Mat4 pij;
    Mat5 mijk;
    float m_step_size;

    MpDca();
	MpDca(std::string mol_type, float pw);
	void init_val();
    void cal_pi();
    void cal_pij();
    void solve_pi();
    void solve_pij();
    void update_eij(float &diff);
    void update_hi(float &diff);
    virtual void calculate_eij();
    virtual float cal_di(int i, int j);
    virtual void set_step(float);
};

} // namespace dca
} // namespace jian


