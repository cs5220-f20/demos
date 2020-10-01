#!/bin/bash
source /etc/profile.d/modules.sh
module load openmpi-4.0.0
mpirun -np 7 $HOME/demos/graphite-hello/mpi-hello.x
