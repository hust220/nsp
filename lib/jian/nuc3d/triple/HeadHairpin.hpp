#include "Module.hpp"

namespace jian {
namespace triple {

class HeadHairpin : public Module {
public:
    HeadHairpin(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace triple
} // namespace jian


