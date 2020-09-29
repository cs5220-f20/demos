CC=gcc-10
CFLAGSA= -fopenmp
CFLAGSB= $(CFLAGSA) -O3 -march=native
CFLAGSC= $(CFLAGSB) -ffast-math

all: ver0a.x ver0b.x ver0c.x ver1c.x ver2c.x ver3c.x

ver-1a.x: kernel-1.c
	$(CC) $(CFLAGSA) $< -o $@

ver0a.x: kernel0.c
	$(CC) $(CFLAGSA) $< -o $@

ver0b.x: kernel0.c
	$(CC) $(CFLAGSB) $< -o $@

ver0c.x: kernel0.c
	$(CC) $(CFLAGSC) $< -o $@

ver1c.x: kernel1.c
	$(CC) $(CFLAGSC) $< -o $@

ver2c.x: kernel2.c
	$(CC) $(CFLAGSC) $< -o $@

ver3c.x: kernel3.c
	$(CC) $(CFLAGSC) $< -o $@
