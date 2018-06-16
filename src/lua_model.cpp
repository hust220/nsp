#include "lua_model.hpp"
#include "pdb.hpp"

namespace jian {

static int lua_pushchain(lua_State *L, Chain *p) {
    Chain ** r = new Chain *;
    *r = p;
    lua_pushlightuserdata(L, r);

    luaL_getmetatable(L, "Chain.meta");
    lua_setmetatable(L, -2);

}

static int model_new(lua_State *L) {
    Model **c = (Model **)lua_newuserdata(L, sizeof(Model *));
    *c = new Model;

    if (lua_gettop(L) == 2) {
        mol_read(**c, lua_tolstring(L, 1, NULL));
    }

    luaL_getmetatable(L, "Model.meta");
    lua_setmetatable(L, -2);

    return 1;
}

static int model_gc (lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, 1);
    if (*c) delete *c;
    return 1;
}

static int model_iter_closure(lua_State *L) {
    Model *c = *(Model **)lua_touserdata(L, lua_upvalueindex(1));
    Int i = lua_tointeger(L, lua_upvalueindex(2));
    if (i < size(*c)) {
        lua_pushchain(L, &(c->at(i)));

        i++;
        lua_pushinteger(L, i);
        lua_replace(L, lua_upvalueindex(2));
    }
    else {
        lua_pushnil(L);
    }
    return 1;
}

static int model_iter(lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, 1);
    lua_pushlightuserdata(L, c);
    lua_pushinteger(L, 0);
    lua_pushcclosure(L, &model_iter_closure, 2);
    return 1;
}

static int model_tostring(lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, -1);
    std::stringstream stream;
    stream << **c;
    lua_pushstring(L, stream.str().c_str());
    return 1;
}

static int model_read(lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, 1);
    auto filename = lua_tolstring(L, 2, NULL);
    mol_read(**c, filename);
    return 1;
}

static int model_getname(lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, -1);
    lua_pushstring(L, (*c)->name.c_str());
    return 1;
}

static int model_setname(lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, -2);
    auto value = lua_tolstring(L, -1, NULL);
    (*c)->name = value;
    return 1;
}

static int model_push(lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, -2);
    Chain **r = (Chain **)lua_touserdata(L, -1);
    (*c)->push_back(**r);
    return 1;
}

static int model_get(lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, 1);
    auto ind = lua_tointeger(L, 2);
    lua_pushchain(L, &((*c)->at(ind-1)));
    return 1;
}

static int model_set(lua_State *L) {
    Model **c = (Model **)lua_touserdata(L, 1);
    auto ind = lua_tointeger(L, 2);
    Chain **r = (Chain **)lua_touserdata(L, 3);
    (*c)->at(ind) = **r;
    return 1;
}

static const struct luaL_Reg modellib_f [] = {
    {"new", model_new},
    {NULL, NULL}
};

static const struct luaL_Reg modellib_m [] = {
    {"__gc", model_gc},
    {"__tostring", model_tostring},
    {"iter", model_iter},
    {"read", model_read},
    {"getname", model_getname},
    {"setname", model_setname},
    {"push", model_push},
    {"get", model_get},
    {"set", model_set},
    {NULL, NULL}
};

int lua_openmodel(lua_State *L) {
    luaL_newmetatable(L, "Model.meta");
    lua_pushstring(L, "__index");
    lua_pushvalue(L, -2); /* pushes the metatable */
    lua_settable(L, -3); /* metatable.__index = metatable */
    luaL_setfuncs(L, modellib_m, 0);

    luaL_newlib(L, modellib_f);
    lua_setglobal(L, "Model");
    return 1;
}

}

