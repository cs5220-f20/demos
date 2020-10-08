MPICC=mpicc

.PHONY: all clean

all: telephone

telephone: telephone.o garble.o
	$(MPICC) -o telephone $^

telephone.o: telephone.c garble.h
garble.o: garble.c garble.h

%.o: %.c
	$(MPICC) -c $<

clean:
	rm -f telephone telephone.o garble.o
