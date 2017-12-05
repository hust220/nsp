#include "dca_plm.hpp"

BEGIN_JN

namespace dca {

    REG_DCA_FAC("plm", PlmDca);

    PlmDca::PlmDca() {}

    PlmDca::PlmDca(S mol_type, float pw) {}

    void PlmDca::calculate_eij() {}

    float PlmDca::cal_di(int i, int j) {return 0;}


}

END_JN

