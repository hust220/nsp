#pragma once

#include <iostream>
#include <map>
#include <functional>
#include <jian/utils/Par.hpp>
#include <jian/utils/file.hpp>
#include <jian/utils/Env.hpp>
#include <jian/utils/log.hpp>
#include <string>
#include <boost/preprocessor.hpp>

namespace jian {

#define REGISTER_NSP_COMPONENT(i)\
    void BOOST_PP_CAT(nsp_, i)(Par par);\
    NSPComponent BOOST_PP_CAT(component_, i)(#i, BOOST_PP_CAT(nsp_, i));\
    void BOOST_PP_CAT(nsp_, i)(Par par)

class NSP {
public:
    std::map<std::string, std::function<void(Par)>> _methods;

    static NSP &instance();

    static void run(int argc, char **argv);

};

class NSPComponent : public NSP {
public:
    template<typename F>
    NSPComponent(const std::string &name, F &&f) {
        NSP::instance()._methods[name] = f;
    }
};

} // namespace jian

