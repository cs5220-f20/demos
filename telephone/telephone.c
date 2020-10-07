#include "garble.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
 * HW2: Telephone
 *
 * Most of the code is already written; you should try to complete the
 * assigment only by completing the lines indicated. Please note that
 * this assignment may be finished with the addition of only two lines of
 * code. However, if you feel it necessary, you may modify more than just
 * the completion areas.
 */

int main(int argc, char **argv)
{
    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    if(argc != 2) {
        printf("Telephone requires a message: \n \
            Please pass in a nonempty string into the command line");
        MPI_Finalize();
        return 0;
    }
    char *buf = argv[1];
    int len = strlen(buf);

    // Get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Seed Random Number Generator, making sure ranks have different seeds
    srand((unsigned)(clock() + world_rank));

    // Start Telephone
    if (world_rank==0) {
        printf("Starting Telephone with %d MPI Ranks...\n", world_size);
        printf("MPI rank 0 starting message: %s \n", buf);
    }

    // If there is just one rank, bail out
    if (world_size == 1) {
        printf("Only one process, refusing to talk to self.\n");
        MPI_Finalize();
        return 0;
    }

    // Send and Receive Loop
    int i;
    for (i = 0; i < world_size; i++) {
        if (world_rank == i) {

            /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ~~~~~~~~~~~~~~Complete: ~~~~~~~~~~~~~~~~

            MPI_Send();

            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

        } else if (world_rank == (i+1) % world_size) {

            /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ~~~~~~~~~~~~~~Complete: ~~~~~~~~~~~~~~~~

            MPI_Recv();

            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

            garble(buf);
            printf("MPI rank %d received message: %s\n", world_rank, buf);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Finalize MPI environment.
    MPI_Finalize();
    return 0;
}
