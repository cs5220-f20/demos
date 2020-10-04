CC=mpicc
CFLAGS=-Wall
MPIRUN=mpirun

.PHONY: run-local run-cluster clean

run-local: laptop.txt

laptop.txt: ping.x
	$(MPIRUN) -n 2 ./ping.x > laptop.txt

run-cluster: ping.x
	sbatch --requeue ping.sub

ping.x: ping.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ping.x laptop.txt 

