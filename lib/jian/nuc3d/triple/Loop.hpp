#include "Module.hpp"

namespace jian {
namespace triple {

class Loop : public Module {
public:
    Loop(const Tuple &, const Tuple &);
    virtual std::string type() const;
};

} // namespace triple
} // namespace jian


