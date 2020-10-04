#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv)
{
    int world_size;
    int world_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    printf("Hello from %d / %d\n", world_rank, world_size);
    MPI_Finalize();
    return 0;
}
