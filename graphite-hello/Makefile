.PHONY: run1 run2 clean

mpi-hello.x: mpi-hello.c
	mpicc -o $@ $<

run1:
	sbatch --requeue mpi_hello_world.sub

run2: mpi-hello.x
	sbatch --requeue mpi_hello_world2.sub

clean:
	rm mpi-hello.x
