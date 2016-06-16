#include "Module.hpp"

namespace jian {
namespace triple {

class TailHairpin : public Module {
public:
    TailHairpin(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace triple
} // namespace jian


