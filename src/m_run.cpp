#include "nsp.hpp"
#include "lua.hpp"
#include <nsp/pdb.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(run) {
    auto g = par.getv("global");
    lua_run(g[1], par);
}

REGISTER_NSP_COMPONENT(cmd) {
    auto g = par.getv("global");
    lua_cmd(g[1], par);
}

END_JN

