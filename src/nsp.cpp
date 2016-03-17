#include "nsp.hpp"
#include <jian/etl/exception.hpp>

int main(int argc, char **argv) {
    try {
        jian::NSP::run(argc, argv);
    } catch (const jian::Error &inf) {
        std::cout << inf.what() << std::endl;
    } catch (const char * inf) {
        std::cout << inf << std::endl;
    } catch (const std::string &s) {
        std::cout << s << std::endl;
    }
}

