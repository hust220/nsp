#include "lua_atom.hpp"
#include <nsp/pdb.hpp>
#include <jian/geom.hpp>

BEGIN_JN

static int atom_new(lua_State *L) {
    Str name = lua_tostring(L, 1);
    Num x = lua_tonumber(L, 2);
    Num y = lua_tonumber(L, 3);
    Num z = lua_tonumber(L, 4);
    Atom **a = (Atom **)lua_newuserdata(L, sizeof(Atom *));
    *a = new Atom(name, x, y, z);

    luaL_getmetatable(L, "Atom.meta");
    lua_setmetatable(L, -2);

    lua_reg(L, a, true);

    return 1; /* new userdatum is already on the stack */
}

static int atom_gc (lua_State *L) {
    lua_gc<Atom>(L);
    return 1;
}

static int atom_tostring(lua_State *L) {
    Atom **a = (Atom **)lua_touserdata(L, 1);
    std::stringstream stream;
    stream << **a;
    lua_pushstring(L, stream.str().c_str());
    return 1;
}

static int atom_dist(lua_State *L) {
    Atom **a1 = (Atom **)lua_touserdata(L, 1);
    Atom **a2 = (Atom **)lua_touserdata(L, 2);
    lua_pushnumber(L, geom::distance(**a1, **a2));
    return 1;
}

static int atom_dist2(lua_State *L) {
    Atom **a1 = (Atom **)lua_touserdata(L, 1);
    Atom **a2 = (Atom **)lua_touserdata(L, 2);
    lua_pushnumber(L, geom::dist2(**a1, **a2));
    return 1;
}

static int atom_ang(lua_State *L) {
    Atom **a1 = (Atom **)lua_touserdata(L, 1);
    Atom **a2 = (Atom **)lua_touserdata(L, 2);
    Atom **a3 = (Atom **)lua_touserdata(L, 3);
    lua_pushnumber(L, geom::angle(**a1, **a2, **a3));
    return 1;
}

static int atom_dih(lua_State *L) {
    Atom **a1 = (Atom **)lua_touserdata(L, 1);
    Atom **a2 = (Atom **)lua_touserdata(L, 2);
    Atom **a3 = (Atom **)lua_touserdata(L, 3);
    Atom **a4 = (Atom **)lua_touserdata(L, 4);
    lua_pushnumber(L, geom::dihedral(**a1, **a2, **a3, **a4));
    return 1;
}

static int atom_get(lua_State *L) {
    Atom **a = (Atom **)lua_touserdata(L, 1);
    if (lua_isinteger(L, 2)) {
        Int ind = lua_tointeger(L, 2);
        lua_pushnumber(L, (*a)->at(ind-1));
    }
    else if (lua_isstring(L, 2)) {
        Str name = lua_tostring(L, 2);
        if (name == "name") {
            lua_pushstring(L, (*a)->name.c_str());
        }
    }
    return 1;
}

static int atom_set(lua_State *L) {
    Atom **a = (Atom **)lua_touserdata(L, 1);
    if (lua_isinteger(L, 2)) {
        Int ind = lua_tointeger(L, 2);
        Num value = lua_tonumber(L, 3);
        (*a)->at(ind-1) = value;
    }
    else if (lua_isstring(L, 2)) {
        Str name = lua_tostring(L, 2);
        Str value = lua_tostring(L, 3);
        if (name == "name") {
            (*a)->name = value;
        }
    }
    return 1;
}

static const struct luaL_Reg atomlib_f [] = {
    {"new", atom_new},
    {"dist", atom_dist},
    {"dist2", atom_dist2},
    {"ang", atom_ang},
    {"dih", atom_dih},
    {NULL, NULL}
};

static const struct luaL_Reg atomlib_m [] = {
    {"__gc", atom_gc},
    {"__tostring", atom_tostring},
    {"__index", atom_get},
    {"__newindex", atom_set},
    {NULL, NULL}
};

int lua_openatom(lua_State *L) {
    luaL_newmetatable(L, "Atom.meta");
//    lua_pushstring(L, "__index");
//    lua_pushvalue(L, -2); /* pushes the metatable */
//    lua_settable(L, -3); /* metatable.__index = metatable */
    luaL_setfuncs(L, atomlib_m, 0);

    luaL_newlib(L, atomlib_f);
//    lua_getglobal(L, "Atom");
//    if (lua_isnil(L, -1)) {
//        lua_pop(L, 1);
//        lua_newtable(L);
//    }
//    luaL_setfuncs(L, atomlib, 0);
    lua_setglobal(L, "Atom");
    return 1;
}

END_JN

