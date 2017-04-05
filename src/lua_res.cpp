#include "lua_res.hpp"
#include <nsp/pdb.hpp>

BEGIN_JN

static int res_new(lua_State *L) {
    Residue **r = (Residue **)lua_newuserdata(L, sizeof(Residue *));
    *r = new Residue;

    luaL_getmetatable(L, "Residue.meta");
    lua_setmetatable(L, -2);

    lua_reg(L, r, true);

    return 1;
}

static int res_gc (lua_State *L) {
    lua_gc<Residue>(L);
    return 1;
}

static int lua_pushatom(lua_State *L, Atom *a) {
    //Atom ** r = (Atom **)std::malloc(sizeof(Atom *));
    Atom ** r = (Atom **)lua_newuserdata(L, sizeof(Atom *));
    //Atom ** r = new Atom *;
    *r = a;
    //lua_pushlightuserdata(L, r);

    luaL_getmetatable(L, "Atom.meta");
    lua_setmetatable(L, -2);

    lua_reg(L, r, false);
}

static int res_iter_closure(lua_State *L) {
    Residue *c = *(Residue **)lua_touserdata(L, lua_upvalueindex(1));
    Int i = lua_tointeger(L, lua_upvalueindex(2));
    if (i < size(*c)) {
        lua_pushatom(L, &(c->at(i)));

        i++;
        lua_pushinteger(L, i);
        lua_replace(L, lua_upvalueindex(2));
    }
    else {
        lua_pushnil(L);
    }
    return 1;
}

static int res_iter(lua_State *L) {
    //Residue **c = (Residue **)lua_touserdata(L, 1);
    lua_pushvalue(L, 1);
    lua_pushinteger(L, 0);
    lua_pushcclosure(L, &res_iter_closure, 2);
    return 1;
}

static int res_tostring(lua_State *L) {
    Residue **a = (Residue **)lua_touserdata(L, -1);
    std::stringstream stream;
    stream << **a;
    lua_pushstring(L, stream.str().c_str());
    return 1;
}

static int res_push(lua_State *L) {
    Residue **r = (Residue **)lua_touserdata(L, -2);
    Atom **a = (Atom **)lua_touserdata(L, -1);
    (*r)->push_back(**a);
    return 1;
}

static int res_get(lua_State *L) {
    Residue **r = (Residue **)lua_touserdata(L, 1);
    if (lua_isinteger(L, 2)) {
      Int ind = lua_tointeger(L, 2);
      lua_pushatom(L, &((*r)->at(ind-1)));
    }
    else if (lua_isstring(L, 2)) {
        Str name = lua_tostring(L, 2);
        if (name == "name") {
            lua_pushstring(L, (*r)->name.c_str());
        }
        else if (name == "len") {
            lua_pushinteger(L, size(**r));
        }
        else {
            lua_pushatom(L, &((**r)[Str(name)]));
        }
    }
    else {
        std::cout << "hi" << std::endl;
    }
    return 1;
}

static int res_set(lua_State *L) {
    Residue **r = (Residue **)lua_touserdata(L, 1);
    Atom **a = (Atom **)lua_touserdata(L, 3);

    if (lua_isinteger(L, 2)) {
        Int ind = lua_tointeger(L, 2);
        (**r)[ind-1] = **a;
    }
    else if (lua_isstring(L, 2)) {
        Str name = lua_tostring(L, 2);
        if (name == "name") {
            (*r)->name = name;
        }
        else {
            (**r)[name] = **a;
        }
    }
    return 1;
}

static const struct luaL_Reg reslib_f [] = {
    {"new", res_new},
    {"iter", res_iter},
    {"push", res_push},
    {NULL, NULL}
};

static const struct luaL_Reg reslib_m [] = {
    {"__gc", res_gc},
    {"__index", res_get},
    {"__newindex", res_set},
    {"__tostring", res_tostring},
    {NULL, NULL}
};

int lua_openres(lua_State *L) {
    luaL_newmetatable(L, "Residue.meta");
    luaL_setfuncs(L, reslib_m, 0);

    luaL_newlib(L, reslib_f);
    lua_setglobal(L, "Residue");
    return 1;
}

END_JN

