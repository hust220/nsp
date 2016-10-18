#include "nsp.hpp"
#include <jian/mcsm/MCpsb.hpp>

namespace jian {
namespace {

class EnPsb : public nuc3d::mc::MCpsb {
    virtual void mc_select() {}
    virtual bool is_selected(const int &i) const {return true;}
    virtual Vec rotating_center() const {return Vec{};}
};

REGISTER_NSP_COMPONENT(en_psb) {
    EnPsb en_psb;
    en_psb.init(par);
    en_psb.total_energy();
}

}
} // namespace jian
















