#include "BuildJunction.h"

int main() {
    try {
        jian::nuc3d::BuildJunction<jian::pdb::RNA> build_junction;
        build_junction("(((())(())))");
    } catch (const char *s) {
        std::cout << s << std::endl;
    }
    return 0;
}

