#include "lua_chain.hpp"
#include "pdb.hpp"

BEGIN_JN

static int lua_pushres(lua_State *L, Residue *p) {
    //Residue ** r = new Residue *;
    Residue ** r = (Residue **)lua_newuserdata(L, sizeof(Residue *));
    *r = p;
    //lua_pushlightuserdata(L, r);

    luaL_getmetatable(L, "Residue.meta");
    lua_setmetatable(L, -2);

    lua_reg(L, r, false);

}

static int chain_new(lua_State *L) {
    Chain **c = (Chain **)lua_newuserdata(L, sizeof(Chain *));
    *c = new Chain;

    if (lua_gettop(L) == 2) {
        chain_read_model(**c, lua_tolstring(L, 1, NULL));
    }

    luaL_getmetatable(L, "Chain.meta");
    lua_setmetatable(L, -2);

    lua_reg(L, c, true);

    return 1;
}

static int chain_gc (lua_State *L) {
    lua_gc<Chain>(L);
    return 1;
}

static int chain_iter_closure(lua_State *L) {
    Chain *c = *(Chain **)lua_touserdata(L, lua_upvalueindex(1));
    Int i = lua_tointeger(L, lua_upvalueindex(2));
    if (i < size(*c)) {
        lua_pushres(L, &(c->at(i)));

        i++;
        lua_pushinteger(L, i);
        lua_replace(L, lua_upvalueindex(2));
    }
    else {
        lua_pushnil(L);
    }
    return 1;
}

static int chain_iter(lua_State *L) {
    //Chain **c = (Chain **)lua_touserdata(L, 1);
    //lua_pushlightuserdata(L, c);
    lua_pushvalue(L, 1);
    lua_pushinteger(L, 0);
    lua_pushcclosure(L, &chain_iter_closure, 2);
    stackDump(L);
    return 1;
}

static int chain_tostring(lua_State *L) {
    Chain **c = (Chain **)lua_touserdata(L, -1);
    std::stringstream stream;
    stream << **c;
    lua_pushstring(L, stream.str().c_str());
    return 1;
}

static int chain_read(lua_State *L) {
    Chain **c = (Chain **)lua_touserdata(L, 1);
    Str filename = lua_tostring(L, 2);
    mol_read(**c, filename);
    return 1;
}

static int chain_push(lua_State *L) {
    Chain **c = (Chain **)lua_touserdata(L, -2);
    Residue **r = (Residue **)lua_touserdata(L, -1);
    (*c)->push_back(**r);
    return 1;
}

static int chain_get(lua_State *L) {
    Chain **c = (Chain **)lua_touserdata(L, 1);
    if (lua_isinteger(L, 2)) {
        Int ind = lua_tointeger(L, 2);
        lua_pushres(L, &((**c)[ind-1]));
    }
    else if (lua_isstring(L, 2)) {
        Str name = lua_tostring(L, 2);
        if (name == "name") lua_pushstring(L, (*c)->name.c_str()); 
        else if (name == "len") lua_pushinteger(L, size(**c));
    }
    return 1;
}

static int chain_set(lua_State *L) {
    Chain **c = (Chain **)lua_touserdata(L, 1);
    if (lua_isinteger(L, 2)) {
        Int ind = lua_tointeger(L, 2);
        Residue **r = (Residue **)lua_touserdata(L, 3);
        (**c)[ind-1] = **r;
    }
    else if (lua_isstring(L, 2)) {
        Str name = lua_tostring(L, 2);
        Str val = lua_tostring(L, 3);
        if (name == "name") {
            (*c)->name = val;
        }
    }
    return 1;
}

static const struct luaL_Reg chainlib_f [] = {
    {"new", chain_new},
    {"iter", chain_iter},
    {"read", chain_read},
    {"push", chain_push},
    {NULL, NULL}
};

static const struct luaL_Reg chainlib_m [] = {
    {"__gc", chain_gc},
    {"__tostring", chain_tostring},
    {"__index", chain_get},
    {"__newindex", chain_set},
    {NULL, NULL}
};

int lua_openchain(lua_State *L) {
    luaL_newmetatable(L, "Chain.meta");
    luaL_setfuncs(L, chainlib_m, 0);

    luaL_newlib(L, chainlib_f);
    lua_setglobal(L, "Chain");
    return 1;
}

END_JN

