CC=gcc
CFLAGS=-O3 -ffast-math -march=native 
LIBS= 

wave_demo.x: wave_demo.c
	$(CC) $(CFLAGS) -o wave_demo.x wave_demo.c $(LIBS)

%.md: %.c
	ldoc wave_demo.c -o wave_demo.md

%.pdf: %.md
	pandoc $< -o $@

%.html: %.md
	pandoc --mathjax --standalone $< -o $@

clean:
	rm -f wave_demo.x
