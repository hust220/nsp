#pragma once

#include <iostream>
#include <map>
#include <functional>
#include <jian/etl/util/Par.hpp>

namespace jian {

#define REGISTER_NSP_COMPONENT(i)\
    void BOOST_PP_CAT(nsp_, i)(const Par &par);\
    NSPComponent BOOST_PP_CAT(component_, i)(#i, BOOST_PP_CAT(nsp_, i));\
    void BOOST_PP_CAT(nsp_, i)(const Par &par)

class NSP {
public:
    std::map<std::string, std::function<void(const Par &)>> _methods;

    static NSP &instance() {
        static NSP nsp;
        return nsp;
    }

    static void run(int argc, char **argv) {
        Par par(argc, argv);
        instance()._methods[par[1]](par);
    }

};

class NSPComponent : public NSP {
public:
    template<typename F>
    NSPComponent(const std::string &name, F &&f) {
        NSP::instance()._methods[name] = f;
    }
};

} // namespace jian

