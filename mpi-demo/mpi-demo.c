#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// Hello world
void hello(int size, int rank)
{
    printf("Hello from %d / %d\n", rank, size);
}

// Send and recv - this deadlocks when the message size gets big
void deadlock_attempt(int size, int rank)
{
    if (size != 2) {
        printf("Should be run with two ranks\n");
	return;
    }

    int num_ints = 1;
    while (1) {
        MPI_Status status;
        int* sendbuf = (int*) malloc(num_ints * sizeof(int));
        int* recvbuf = (int*) malloc(num_ints * sizeof(int));
	assert(sendbuf && recvbuf);
        MPI_Send(sendbuf, num_ints, MPI_INT, 1-rank, 0, MPI_COMM_WORLD);
        MPI_Recv(recvbuf, num_ints, MPI_INT, 1-rank, 0, MPI_COMM_WORLD, &status);
        printf("%d: Received from %d (%d)\n", rank, status.MPI_SOURCE, num_ints);
        free(recvbuf);
        free(sendbuf);
	num_ints *= 10;
    }
}

// Sendrecv - this does not deadlock, even for big message sizes
void deadlock_nonattempt(int size, int rank)
{
    if (size != 2) {
        printf("Should be run with two ranks\n");
	return;
    }

    int num_ints = 1;
    for (int k = 0; k < 7; ++k) {
        MPI_Status status;
        int* sendbuf = (int*) malloc(num_ints * sizeof(int));
        int* recvbuf = (int*) malloc(num_ints * sizeof(int));
	assert(sendbuf && recvbuf);
        MPI_Sendrecv(sendbuf, num_ints, MPI_INT, 1-rank, 0, 
                     recvbuf, num_ints, MPI_INT, 1-rank, 0, 
		     MPI_COMM_WORLD, &status);
        printf("%d: Received from %d (%d)\n", rank, status.MPI_SOURCE, num_ints);
        free(recvbuf);
        free(sendbuf);
	num_ints *= 10;
    }
}

// Scatter data and use nonblocking communication for information exchange
void compute_interacts(int size, int rank)
{
    // Set up a random array of points (rank 0)
    double* src_array = NULL;
    int ntotal = 5;
    if (rank == 0) {
        src_array = (double*) malloc(ntotal * sizeof(double));
	printf("Initialize array:");
	for (int i = 0; i < ntotal; ++i) {
	    src_array[i] = drand48();
	    printf(" %g", src_array[i]);
	}
	printf("\n");
    }

    // Set up counts, offsets, and local arrays
    int offsets[size+1];
    int counts[size];
    for (int i = 0; i <= size; ++i) 
        offsets[i] = (i * ntotal)/size;
    for (int i = 0; i < size; ++i)
	counts[i] = offsets[i+1]-offsets[i];
    double* loc_array = (double*) malloc(counts[rank] * sizeof(double));

    // Distribute points to all ranks (scatter)
    MPI_Scatterv(src_array, counts, offsets, MPI_DOUBLE,
		 loc_array, counts[rank], MPI_DOUBLE,
		 0, MPI_COMM_WORLD);

    // Print what we got
    printf("On rank %d (%d-%d): ", rank, offsets[rank], offsets[rank+1]);
    for (int i = 0; i < counts[rank]; ++i)
        printf(" %g", loc_array[i]);
    printf("\n");

    // Compute interactions between all points (n^2), passing data cyclically (TODO)
  
    // Clean up 
    free(loc_array); 
    if (rank == 0)
        free(src_array);
}

// Broadcast
// Reduce
// Scatter/gather
// Scan

int main(int argc, char** argv)
{
    int world_size;
    int world_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Hello world exercise
    hello(world_size, world_rank);

    // Only run this if you want to deadlock! 
    // deadlock_attempt(world_size, world_rank);

    // Or not
    // deadlock_nonattempt(world_size, world_rank);

    // Compute pairwise interactions
    compute_interacts(world_size, world_rank);

    MPI_Finalize();
    return 0;
}
