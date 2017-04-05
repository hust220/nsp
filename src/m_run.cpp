#include "nsp.hpp"
#include "lua.hpp"
#include "lua_atom.hpp"
#include "lua_res.hpp"
#include "lua_chain.hpp"
#include "lua_model.hpp"
#include "lua_mol.hpp"
#include <nsp/pdb.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(run) {
    auto g = par.getv("global");
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);// open lua libraries
    lua_openatom(L);
    lua_openres(L);
    lua_openchain(L);
    lua_openmodel(L);
    lua_openmol(L);
    luaL_dofile(L, g[1].c_str());// run script
    lua_close(L);// close lua
}

REGISTER_NSP_COMPONENT(cmd) {
    auto g = par.getv("global");
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);// open lua libraries
    lua_openatom(L);
    lua_openres(L);
    lua_openchain(L);
    lua_openmodel(L);
    lua_openmol(L);
    luaL_dostring(L, g[1].c_str());// run script
    lua_close(L);// close lua
}

END_JN

