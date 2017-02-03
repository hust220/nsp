#include "nsp.hpp"

BEGIN_JN

std::map<std::string, std::function<void(Par)>> handles
{
    {
        "ddna", [](Par par){
            std::cout << "ddna" << std::endl;
            // 1. assemble
            // Assemble ass(par);
            // ass.predict();
            // 2. if has no templates; then do optimization
        }
    }, 
    {
        "tdna", [](Par par){
            std::cout << "tdna" << std::endl;
            // THMC thmc(par);
            // 1. construct inittial structure
            // thmc.build_initial_scaffold();
            // 2. optimization
        }
    },
    {
        "qdna", [](Par par){
            std::cout << "qdna" << std::endl;
            // 1. construct inittial structure
            // 2. optimization
        }
    }
};

REGISTER_NSP_COMPONENT(dtsp) {
    auto global = par.getv("global");
    Str mode = global[1];
    to_lower(mode);
    handles[mode](par);
}

END_JN

