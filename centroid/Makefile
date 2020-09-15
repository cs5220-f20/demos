CC=gcc
CFLAGS=-fopenmp -std=gnu99 -O2

centroid.x: centroid.c
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f centroid.x
