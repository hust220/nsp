#pragma once

#include <iostream>
#include <map>
#include <functional>
#include <jian/utils/Par.hpp>
#include <jian/utils/file.hpp>
#include <jian/utils/Env.hpp>
#include <jian/utils/log.hpp>
#include <string>
#include <jian/pp.hpp>
//#include <boost/preprocessor.hpp>

namespace jian {

#define REGISTER_NSP_COMPONENT(i)\
    void PP_CAT(nsp_, i)(Par par);\
    NSPComponent  PP_CAT(component_, i)(PP_STRING(i), PP_CAT(nsp_, i));\
    void PP_CAT(nsp_, i)(Par par)

class NSP {
public:
    std::map<str_t, std::function<void(Par)>> _methods;
	int m_argc;
	char **m_argv;

    static NSP &instance();

    static void run(int argc, char **argv);

};

class NSPComponent : public NSP {
public:
    template<typename F>
    NSPComponent(const str_t &name, F &&f) {
        NSP::instance()._methods[name] = f;
    }
};

} // namespace jian

