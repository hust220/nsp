#include "Constraints.h"

int main(int argc, char **argv) {
    jian::dg::Constraints constraints;
    for (int i = 0; i < 3; i++) for (int j = i + 1; j < 5; j++) constraints(i, j) = {double(i + j), 0.0};
    for (auto pair : constraints._len_terms) {
        std::cout << pair.first[0] << '-' << pair.first[1] << ' ' << pair.second.avg << std::endl;
    }
    return 0;
}

