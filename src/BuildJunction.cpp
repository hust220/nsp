#include <jian/nuc3d/BuildJunction.h>
#include <jian/etl.h>

int main(int argc, char **argv) {
    try {
        jian::Par par(argc, argv);
        int n = boost::lexical_cast<int>(par["num"][0]);
        jian::nuc3d::BuildJunction<jian::pdb::RNA> build_junction;
        for (int i = 0; i < n; i++) {
            write_pdb(build_junction(par["ss"][0]), par["name"][0]+'-'+boost::lexical_cast<std::string>(i+1)+".pdb");
        }
    } catch (const char *s) {
        std::cout << s << std::endl;
    }
    return 0;
}


