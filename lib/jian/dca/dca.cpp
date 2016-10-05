#include "dca.hpp"
#include "MfDca.hpp"

namespace jian {
namespace dca {

void analyze(std::string fa_file, std::string out_file, int n) {
    Dca *dca = new MfDca();
    dca->run(fa_file, out_file, n);
    delete dca;
}

} // namespace analyze
} // namespace jian


