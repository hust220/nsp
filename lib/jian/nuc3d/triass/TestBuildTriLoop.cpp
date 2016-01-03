#include "BuildTriLoop.h"
#include "../../pdb/DNA.h"

int main(int argc, char **argv) {
    jian::nuc3d::triass::BuildTriLoop<jian::DNA> build_tri_loop;
    build_tri_loop(std::vector<int>{2, 2, 2}, 0, 1);
    return 0;
}
