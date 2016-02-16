#include "Cif.h"
#include "../etl.h"

int main(int argc, char **argv) {
    jian::Par par(argc, argv);
    auto file_name = par["global"][0];
    jian::Cif cif(file_name);
    return 0;
}


