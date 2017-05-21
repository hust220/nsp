#pragma once

#include <iostream>
#include <map>
#include <functional>
#include "par.hpp"
#include "file.hpp"
#include "env.hpp"
#include "log.hpp"
#include <string>
#include "pp.hpp"

BEGIN_JN

#define REGISTER_NSP_COMPONENT(i)\
    void PP_CAT(nsp_, i)(JN_ Par par);\
    NSPComponent PP_CAT(component_, i)(PP_STRING(i), PP_CAT(nsp_, i));\
    void PP_CAT(nsp_, i)(JN_ Par par)

#define ALIAS_NSP_COMPONENT(i, j)\
	NSPComponent PP_CAT(component_, j)(PP_STRING(j), PP_CAT(nsp_, i))

class NSP {
public:
    Ms<Fn<void(Par)>> _methods;
	int m_argc;
	char **m_argv;

    static NSP &instance();

    static void run(int argc, char **argv);

};

class NSPComponent : public NSP {
public:
    template<typename F>
    NSPComponent(const Str &name, F &&f) {
        NSP::instance()._methods[name] = f;
    }
};

END_JN

