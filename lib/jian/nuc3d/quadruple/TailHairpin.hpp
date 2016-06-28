#include "Module.hpp"

namespace jian {
namespace nuc3d {
namespace quadruple {

class TailHairpin : public Module {
public:
    TailHairpin(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace quadruple
}
} // namespace jian


