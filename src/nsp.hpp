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
//    std::map<std::string, std::function<>> _synopses;
//    std::map<std::string, std::string> _details;

    static NSP &instance() {
        static NSP nsp;
        return nsp;
    }

    static void run(int argc, char **argv) {
        Par par(argc, argv);
        if (par.has("log_level")) {
            set_log_level(std::stoi(par["log_level"][0]));
        }
        auto &m = instance()._methods;
        std::string path = Env::lib() + "/RNA/pars/src/";
        if (m.find(par[1]) == m.end()) {
            std::string name = path + "nsp.md";
            EACH_LINE(name.c_str(), std::cout << L << std::endl;);
            for (auto && pair : m) {std::cout << pair.first << ' ';}
            std::cout << std::endl;
        } else if (par.has("help") || par.has("h") || par.has("-help")) {
            std::string name = path + par[1] + ".md";
            EACH_LINE(name.c_str(), std::cout << L << std::endl;);
        } else {
            m[par[1]](par);
        }
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

