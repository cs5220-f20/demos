//ldoc on
/**
 * # Lua interface
 * 
 * It is often convenient to be able to change the parameters for
 * a simulation dynamically, without having to recompile.  Moreover,
 * fiddling with the parameters may be the type of thing that's most
 * easily done in a higher level language than C or Fortran.  This is
 * an example of one way to provide such functionality by embedding
 * an interpreter in the code.  In this case, we use the Lua language,
 * which is widely used as an embedded scripting and configuration
 * language in computer games (also used in TeX, Pandoc, and a surprising
 * number of other tools), but originated as a configuration language
 * for finite element simulation codes. 
 * The [Lua reference manual](http://www.lua.org/manual/5.3) is
 * pretty good for figuring out how to get around the language and the
 * C interfaces, as is [the book](https://www.lua.org/pil/).
 * 
 * You don't strictly *need* to understand anything in this file. It is
 * provided purely in the hopes that you will find it an interesting
 * (and potentially useful) programming pattern.
 * 
 */
//ldoc off

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#include "wave1d.h"

//ldoc on
/**
 * ## Wave stepper bindings in Lua
 *
 * For constructing the simulation, we use named (rather than positional)
 * arguments.  In Lua, these are conventionally passed through a table
 * argument.  The arguments are:
 * 
 * - `n` (1000): Number of mesh points, including boundary
 * - `c` (1.0): Wave speed
 * - `dt`: Time step (default so that the CFL number is 0.8)
 *
 */
static
int wave1dL_alloc(lua_State* L)
{
    if (lua_gettop(L) != 1)
        return luaL_error(L, "expecting exactly one table argument");
    luaL_checktype(L, 1, LUA_TTABLE);
    lua_getfield(L, 1, "n");
    lua_getfield(L, 1, "c");
    lua_getfield(L, 1, "dt");
    int n    = luaL_optinteger(L, 2, 1000);
    float c  = luaL_optnumber(L, 3, 1.0);
    float dt = luaL_optnumber(L, 4, 0.8/c/(n-1));
    wave1d_t** psim = lua_newuserdata(L, sizeof(wave1d_t*));
    *psim = wave1d_alloc(n, c, dt);
    luaL_getmetatable(L, "wave1d.sim");
    lua_setmetatable(L, -2);
    return 1;
}

/**
 * The free function maps to the garbage collection meta-method, but we
 * also allow it to be called explicitly.
 * 
 */
static
int wave1dL_free(lua_State* L)
{
    if (lua_gettop(L) != 1)
        return luaL_error(L, "exactly one argument required");
    wave1d_t** psim = (wave1d_t**) luaL_checkudata(L, 1, "wave1d.sim");
    if (*psim)
        wave1d_free(*psim);
    *psim = NULL;
    return 0;
}

/**
 * We use the wave1dL_check wrapper to ensure that accesses to a deleted
 * sim object result in an error.
 *
 */
static
wave1d_t* wave1dL_check(lua_State* L, int idx)
{
    wave1d_t* sim = *(wave1d_t**) luaL_checkudata(L, idx, "wave1d.sim");
    luaL_argcheck(L, sim != NULL, idx, "Attempt to use deleted sim");
    return sim;
}

/**
 * We define a getter meta-method (via __index)
 * to get stuff out of the wave struct (current step or properties).
 *
 */
static
int wave1dL_getter(lua_State* L)
{
    wave1d_t* sim = wave1dL_check(L, 1);

    // If integer index, read from current frame
    if (lua_isnumber(L, 2)) {
        int idx = luaL_checknumber(L, 2);
        luaL_argcheck(L, 0 <= idx && idx < sim->n, 2, "Index out of bounds");
        lua_pushnumber(L, wave1d_frame(sim, sim->tidx)[idx]);
        return 1;
    }

    // If string index, get a property
    const char* field = luaL_checkstring(L, 2);
    if (strcmp(field, "n") == 0)
        lua_pushinteger(L, sim->n);
    else if (strcmp(field, "C") == 0)
        lua_pushnumber(L, sqrtf(sim->C2));
    else if (strcmp(field, "c") == 0)
        lua_pushnumber(L, sim->c);
    else if (strcmp(field, "dt") == 0)
        lua_pushnumber(L, sim->dt);
    else if (strcmp(field, "tidx") == 0)
        lua_pushinteger(L, sim->tidx);
    else if (strcmp(field, "t") == 0)
        lua_pushnumber(L, sim->tidx * sim->dt);
    else if (strcmp(field, "batch") == 0)
        lua_pushnumber(L, wave1d_batch());
    else { // Otherwise, check metatable
        luaL_getmetatable(L, "wave1d.sim");
        lua_pushstring(L, field);
        lua_gettable(L, -2);
    }
    return 1;
}

/**
 * The time stepper interface takes the simulation and nsteps and advances
 * time.
 *
 */
static
int wave1dL_steps(lua_State* L)
{
    wave1d_t* sim = wave1dL_check(L, 1);
    int nsteps = luaL_optinteger(L, 2, 1);
    luaL_argcheck(L, nsteps > 0, 2, "must advance at least one step");
    lua_pushinteger(L, wave1d_steps(sim, nsteps));
    return 1;
}

/**
 * The frame set routine calls a function at each (interior) mesh point
 * in order to get a function value.
 *
 */
static
int wave1dL_set_frame(lua_State* L)
{
    wave1d_t* sim = wave1dL_check(L, 1);
    luaL_checktype(L, 2, LUA_TFUNCTION);
    int tidx = luaL_optinteger(L, 3, sim->tidx);
    float* frame = wave1d_frame(sim, tidx);
    int n = sim->n;
    for (int i = 1; i < n-1; ++i) {
        float x = (1.0*i)/(n-1);
        lua_pushvalue(L, 2);
        lua_pushnumber(L, x);
        lua_call(L, 1, 1);
        frame[i] = lua_tonumber(L,-1);
        lua_pop(L,1);
    }
    return 0;
}

/**
 * The file dump routine writes everything to a file (default name
 * of `frame.txt`).
 * 
 */
static
int wave1dL_dump(lua_State* L)
{
    wave1d_t* sim = wave1dL_check(L, 1);
    const char* fname = luaL_optstring(L, 2, "frame.txt");
    wave1d_dump(sim, fname);
    return 0;
}

// Package global functions
static const struct luaL_Reg wave1dL_gfun[] = {
    {"new",       wave1dL_alloc},
    {NULL, NULL}
};

// Methods (in meta-metatable)
static const struct luaL_Reg wave1dL_mfun[] = {
    {"free",      wave1dL_free},
    {"get",       wave1dL_getter},
    {"step",      wave1dL_steps},
    {"set_frame", wave1dL_set_frame},
    {"dump",      wave1dL_dump},
    {"__index",   wave1dL_getter},
    {"__gc",      wave1dL_free},
    {NULL, NULL}};

int luaopen_wave1d(lua_State* L)
{
    // Set up metatable
    luaL_newmetatable(L, "wave1d.sim");
    luaL_setfuncs(L, wave1dL_mfun, 0);
    lua_pop(L, 1);

    // Set up global method names
    luaL_newlib(L, wave1dL_gfun);
    lua_setglobal(L, "wave1d");    
    return 0;
}
