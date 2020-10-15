LUA_CFLAGS=`pkg-config --cflags lua-5.3`
LUA_LIBS=`pkg-config --libs lua-5.3`

CC=gcc
CFLAGS=-O3 -ffast-math -march=native $(LUA_CFLAGS)
LIBS= $(LUA_LIBS)

wave1d.x: driver.o wave1d.o wave1d_lua.o
	$(CC) $^ -o $@ $(LIBS)

wave_lua.md: wave1d.h wave1d.c wave1d_lua.c
	ldoc wave1d.h wave1d.c wave1d_lua.c -o wave_demo.md

%.pdf: %.md
	pandoc $< -o $@

clean:
	rm -f *.o wave1d.x
