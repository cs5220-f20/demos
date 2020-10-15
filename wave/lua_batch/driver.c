#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#include "wave1d.h"
#include "wave1d_lua.h"

#if defined(_OPENMP)
#include <omp.h>
#else
#include <time.h>
inline double omp_get_wtime() { return (double) clock() / CLOCKS_PER_SEC; }
#endif

static
int lua_wtime(lua_State* L)
{
    lua_pushnumber(L, omp_get_wtime());
    return 1;
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        fprintf(stderr, "Usage: wave1d params.lua\n");
        return -1;
    }
    
    int status = 0;
    lua_State* L = luaL_newstate();
    luaL_openlibs(L);
    luaopen_wave1d(L);
    lua_pushcfunction(L, lua_wtime);
    lua_setglobal(L, "wtime");
    if (luaL_dofile(L, argv[1])) {
        fprintf(stderr, "%s\n", lua_tostring(L,-1));
        status = -1;
    }
    lua_close(L);
    return status;
}
