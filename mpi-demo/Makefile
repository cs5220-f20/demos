.PHONY: run clean
MPICC=mpicc

mpi-demo.x: mpi-demo.c
	$(MPICC) -o $@ $<

run:
	sbatch --requeue mpi-demo.sub

clean:
	rm -f mpi-demo.x
