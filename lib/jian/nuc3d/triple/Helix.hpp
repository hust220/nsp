#include "Module.hpp"

namespace jian {
namespace triple {

class Helix : public Module {
public:
    Helix(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace triple
} // namespace jian


