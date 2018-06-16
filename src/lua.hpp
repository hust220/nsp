#pragma once

#include "jian.hpp"
#include "par.hpp"

extern "C"  
{  
#include "lua.h"  
#include "lualib.h"  
#include "lauxlib.h"  
};

namespace jian {

void stackDump (lua_State *L);

template<typename T>
inline void lua_reg(lua_State *L, T && a, Bool b) {
    lua_pushlightuserdata(L, a);
    lua_pushboolean(L, b);
    lua_settable(L, LUA_REGISTRYINDEX);
}

template<typename T>
inline void lua_gc(lua_State *L, Str msg = "") {
    T **r = (T **)lua_touserdata(L, 1);
    lua_pushlightuserdata(L, r);
    lua_gettable(L, LUA_REGISTRYINDEX);
    if (lua_toboolean(L, -1)) {
        if (msg != "") std::cout << msg << std::endl;
        if (*r) delete *r;
    }
}

void lua_run(Str filename, const Par &par);

void lua_cmd(Str cmd, const Par &par);

}

