#include "Module.hpp"

namespace jian {
namespace nuc3d {
namespace quadruple {

class Loop : public Module {
public:
    Loop(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace quadruple
}
} // namespace jian


