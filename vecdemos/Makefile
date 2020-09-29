CC=gcc

assembly: simple0.s simple1.s simple2.s simple3.s

%.s: %.c
	$(CC) -O3 -fopt-info-vec -fopt-info-vec-missed -march=native -S $<
