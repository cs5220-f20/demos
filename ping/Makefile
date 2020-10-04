CC=mpicc
CFLAGS=-Wall
MPIRUN=mpirun

.PHONY: run-local run-cluster svgs clean

run-local: laptop.txt

laptop.txt: ping.x
	$(MPIRUN) -n 2 ./ping.x > laptop.txt

run-cluster: ping.x
	sbatch --requeue ping.sub

ping.x: ping.c
	$(CC) $(CFLAGS) $< -o $@

svgs: laptop.svg graphite-01.svg graphite-11.svg graphite-23.svg

%.svg: %.txt
	python plotter.py $*

clean:
	rm -f ping.x laptop.txt 

