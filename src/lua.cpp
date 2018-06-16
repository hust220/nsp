#include <cstdio>
#include "lua.hpp"
#include "lua_atom.hpp"
#include "lua_res.hpp"
#include "lua_chain.hpp"
#include "lua_model.hpp"
#include "lua_mol.hpp"
#include "lua_par.hpp"

namespace jian {

void stackDump (lua_State *L) {
    int i;
    int top = lua_gettop(L);
    for (i = 1; i <= top; i++) { /* repeat for each level */
        int t = lua_type(L, i);
        switch (t) {
            case LUA_TSTRING: /* strings */
                printf("`%s'", lua_tostring(L, i));
                break;
            case LUA_TBOOLEAN: /* booleans */
                printf(lua_toboolean(L, i) ? "true" : "false");
                break;
            case LUA_TNUMBER: /* numbers */
                printf("%g", lua_tonumber(L, i));
                break;
            default: /* other values */
                printf("%s", lua_typename(L, t));
                break;
        }
        printf(" "); /* put a separator */
    }
    printf("\n"); /* end the listing */
}

void lua_run(Str filename, const Par &par) {
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);// open lua libraries
    lua_openatom(L);
    lua_openres(L);
    lua_openchain(L);
    lua_openmodel(L);
    lua_openmol(L);
    lua_openpar(L, par);
    luaL_dofile(L, filename.c_str());// run script
    lua_close(L);// close lua
}

void lua_cmd(Str cmd, const Par &par) {
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);// open lua libraries
    lua_openatom(L);
    lua_openres(L);
    lua_openchain(L);
    lua_openmodel(L);
    lua_openmol(L);
    lua_openpar(L, par);
    luaL_dostring(L, cmd.c_str());// run script
    lua_close(L);// close lua
}

}

