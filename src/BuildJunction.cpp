#include <jian/nuc3d/BuildLoop.h>
#include <jian/pdb.h>

int main(int argc, char **argv) {
    try {
        jian::Par par(argc, argv);
        int n = boost::lexical_cast<int>(par["num"][0]);
        jian::nuc3d::BuildLoop build_loop;
        for (int i = 0; i < n; i++) {
            write_pdb(build_loop(par["ss"][0]), par["name"][0]+'-'+boost::lexical_cast<std::string>(i+1)+".pdb");
        }
    } catch (const char *s) {
        std::cout << s << std::endl;
    }
    return 0;
}


