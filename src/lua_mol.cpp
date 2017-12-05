#include "lua_mol.hpp"
#include "pdb.hpp"

BEGIN_JN

static int lua_pushmodel(lua_State *L, Model *p) {
    Model ** r = new Model *;
    *r = p;
    lua_pushlightuserdata(L, r);

    luaL_getmetatable(L, "Model.meta");
    lua_setmetatable(L, -2);

}

static int mol_new(lua_State *L) {
    Molecule **c = (Molecule **)lua_newuserdata(L, sizeof(Molecule *));
    *c = new Molecule;

    if (lua_gettop(L) == 2) {
        mol_read(**c, lua_tolstring(L, 1, NULL));
    }

    luaL_getmetatable(L, "Molecule.meta");
    lua_setmetatable(L, -2);

    return 1;
}

static int mol_gc (lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, 1);
    if (*c) delete *c;
    return 1;
}

static int mol_iter_closure(lua_State *L) {
    Molecule *c = *(Molecule **)lua_touserdata(L, lua_upvalueindex(1));
    Int i = lua_tointeger(L, lua_upvalueindex(2));
    if (i < size(*c)) {
        lua_pushmodel(L, &(c->at(i)));

        i++;
        lua_pushinteger(L, i);
        lua_replace(L, lua_upvalueindex(2));
    }
    else {
        lua_pushnil(L);
    }
    return 1;
}

static int mol_iter(lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, 1);
    lua_pushlightuserdata(L, c);
    lua_pushinteger(L, 0);
    lua_pushcclosure(L, &mol_iter_closure, 2);
    return 1;
}

static int mol_tostring(lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, -1);
    std::stringstream stream;
    stream << **c;
    lua_pushstring(L, stream.str().c_str());
    return 1;
}

static int mol_read(lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, 1);
    auto filename = lua_tolstring(L, 2, NULL);
    mol_read(**c, filename);
    return 1;
}

static int mol_getname(lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, -1);
    lua_pushstring(L, (*c)->name.c_str());
    return 1;
}

static int mol_setname(lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, -2);
    auto value = lua_tolstring(L, -1, NULL);
    (*c)->name = value;
    return 1;
}

static int mol_push(lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, -2);
    Model **r = (Model **)lua_touserdata(L, -1);
    (*c)->push_back(**r);
    return 1;
}

static int mol_get(lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, 1);
    auto ind = lua_tointeger(L, 2);
    lua_pushmodel(L, &((*c)->at(ind-1)));
    return 1;
}

static int mol_set(lua_State *L) {
    Molecule **c = (Molecule **)lua_touserdata(L, 1);
    auto ind = lua_tointeger(L, 2);
    Model **r = (Model **)lua_touserdata(L, 3);
    (*c)->at(ind) = **r;
    return 1;
}

static const struct luaL_Reg mollib_f [] = {
    {"new", mol_new},
    {NULL, NULL}
};

static const struct luaL_Reg mollib_m [] = {
    {"__gc", mol_gc},
    {"__tostring", mol_tostring},
    {"iter", mol_iter},
    {"read", mol_read},
    {"getname", mol_getname},
    {"setname", mol_setname},
    {"push", mol_push},
    {"get", mol_get},
    {"set", mol_set},
    {NULL, NULL}
};

int lua_openmol(lua_State *L) {
    luaL_newmetatable(L, "Molecule.meta");
    lua_pushstring(L, "__index");
    lua_pushvalue(L, -2); /* pushes the metatable */
    lua_settable(L, -3); /* metatable.__index = metatable */
    luaL_setfuncs(L, mollib_m, 0);

    luaL_newlib(L, mollib_f);
    lua_setglobal(L, "Molecule");
    return 1;
}

END_JN

