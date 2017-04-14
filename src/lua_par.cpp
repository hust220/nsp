#include "lua_par.hpp"

BEGIN_JN

static const char * par_name = "par";

static int par_new(lua_State *L) {
    //std::cout << "function: " << __FUNCTION__ << std::endl;
    lua_pushlightuserdata(L, &par_name);
    lua_gettable(L, LUA_REGISTRYINDEX);
    //stackDump(L);
    return 1;
}

static int par_gc(lua_State *L) {
    //std::cout << "function: " << __FUNCTION__ << std::endl;
    Par **p = (Par **)lua_touserdata(L, 1);
    if (*p) delete *p;
    return 1;
}

static int par_get(lua_State *L) {
    //std::cout << "function: " << __FUNCTION__ << std::endl;
    //stackDump(L);
    Par **p = (Par **)lua_touserdata(L, lua_upvalueindex(1));
    Int n = lua_gettop(L);
    for (Int i = 1; i <= n; i++) {
        Str name = lua_tostring(L, i);
        if ((*p)->has(name)) {
            lua_pushstring(L, (*p)->get(name).c_str());
            return 1;
        }
    }
    lua_pushnil(L);
    return 1;
}

static int par_getv(lua_State *L) {
    //std::cout << "function: " << __FUNCTION__ << std::endl;
    Par **p = (Par **)lua_touserdata(L, lua_upvalueindex(1));
    lua_newtable(L);
    Int n = lua_gettop(L);
    for (Int i = 1; i <= n; i++) {
        Str name = lua_tostring(L, i);
        if ((*p)->has(name)) {
            Int ind = 0;
            for (auto && s : (*p)->getv(name)) {
                lua_pushstring(L, s.c_str());
                lua_rawseti(L, -2, ++ind);
            }
            return 1;
        }
    }
    lua_pushnil(L);
    return 1;
}

static int par_index(lua_State *L) {
    //std::cout << "function: " << __FUNCTION__ << std::endl;
    Par **p = (Par **)lua_touserdata(L, 1);
    if (lua_isinteger(L, 2)) {
        Int n = lua_tointeger(L, 2);
        lua_pushstring(L, (**p)[n-1].c_str());
    }
    else if (lua_isstring(L, 2)) {
        Str name = lua_tostring(L, 2);
        if (name == "get") {
            //stackDump(L);
            lua_pushvalue(L, 1);
            lua_pushcclosure(L, &par_get, 1);
        }
        else if (name == "getv") {
            lua_pushvalue(L, 1);
            lua_pushcclosure(L, &par_getv, 1);
        }
    }
    return 1;
}

static const struct luaL_Reg parlib_f [] = {
    {"new", par_new},
    {NULL, NULL}
};

static const struct luaL_Reg parlib_m [] = {
    {"__index", par_index},
    {"__gc", par_gc},
    {NULL, NULL}
};

void lua_openpar(lua_State *L, const Par &par) {
    lua_pushlightuserdata(L, &par_name);
    Par **p = (Par **)lua_newuserdata(L, sizeof(Par *));
    *p = new Par(par);
    //stackDump(L);

    luaL_newmetatable(L, "Par.meta");
    //stackDump(L);
    luaL_setfuncs(L, parlib_m, 0);
    //stackDump(L);

    lua_setmetatable(L, -2);
    //stackDump(L);

    luaL_newlib(L, parlib_f);
    lua_setglobal(L, "Par");

    lua_settable(L, LUA_REGISTRYINDEX);

}

END_JN

