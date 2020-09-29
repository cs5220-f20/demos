CC=gcc-10

all: ver0a.x ver0b.x ver0c.x ver1c.x ver2c.x ver3c.x

ver-1a.x: kernel-1.c
	$(CC) -fopenmp $< -o $@

ver0a.x: kernel0.c
	$(CC) -fopenmp $< -o $@

ver0b.x: kernel0.c
	$(CC) -fopenmp -O3 -march=native $< -o $@

ver0c.x: kernel0.c
	$(CC) -fopenmp -O3 -march=native -ffast-math $< -o $@

ver1c.x: kernel1.c
	$(CC) -fopenmp -O3 -march=native -ffast-math $< -o $@

ver2c.x: kernel2.c
	$(CC) -fopenmp -O3 -march=native -ffast-math $< -o $@

ver3c.x: kernel3.c
	$(CC) -fopenmp -O3 -march=native -ffast-math $< -o $@
