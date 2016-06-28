#include "Module.hpp"

namespace jian {
namespace nuc3d {
namespace quadruple {

class HeadHairpin : public Module {
public:
    HeadHairpin(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace quadruple
}
} // namespace jian


