#include "Module.hpp"

namespace jian {
namespace nuc3d {
namespace quadruple {

class Helix : public Module {
public:
    Helix(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace quadruple
}
} // namespace jian


