#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/**
 * ## Hello world
 *
 * This one is relatively straightforward: each rank says hello by
 * printing to the standard output.  Note that the output lines can
 * appear in whatever order they want!
 *
 */
void hello(int size, int rank)
{
    printf("Hello from %d / %d\n", rank, size);
}


/**
 * ## Send, recv, and sendrecv
 *
 * The `MPI_Send` operation is blocking, which means that under some
 * circumstances (like those we described in class), we can get into a
 * deadlock if there is a cycle of processors trying to send to each
 * other who mutually can't make it to their receive operations until
 * the sends complete.  *But* whether this actually results in deadlock
 * or not depends on the size of the messages and the amount of buffering
 * that the MPI implementation (and operating system) provides!  Here
 * we test this out, running with incrementally larger messages until
 * we actually see a deadlock condition.
 *
 */
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
        MPI_Send(sendbuf, num_ints, MPI_INT, 1-rank,
                 0, MPI_COMM_WORLD);
        MPI_Recv(recvbuf, num_ints, MPI_INT, 1-rank,
                 0, MPI_COMM_WORLD, &status);
        printf("%d: Received from %d (%d)\n",
               rank, status.MPI_SOURCE, num_ints);
        free(recvbuf);
        free(sendbuf);
        num_ints *= 10;
    }
}


/**
 * The pairing of an overlapping send and receive is pretty common,
 * enough to deserve special support with `MPI_Sendrecv`.  Though
 * the overall send-receive pair blocks, the sub-operations (the
 * individual sends and recieves) run concurrently, so that we won't
 * deadlock even if there is a communication cycle.  We test this
 * out here.
 *
 */
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


/**
 * ## Computing all-to-all interactions
 *
 * The typical all-to-all interaction looks like
 * $$
 *   f(x_i) = \sum_{j} k(x_i, x_j) c_j
 * $$
 * where $k$ is an interaction kernel.  For some kernels, there is
 * a faster way to compute this than explicitly evaluating all the
 * sums, but that is not our interest here.  We want to show off
 * some MPI calls by doing things in a brute-force way.
 *
 * The "pass it around" protocol discussed in the slides manages
 * memory scalability by saying that each rank is assigned "ownership"
 * of disjoint pieces of the array of points $x_i$.  A rank is responsible
 * for accumulating $f(x_i)$ for the points that it owns; but to do so,
 * it needs to see all the other points.  We manage this with a
 * "pass the dish" protocol: we keep a set of particle data that we
 * compute interactions against at each phase, and rotate those interactions
 * around the ring.
 *
 * In class, we wrote this as one long function.  I'm less against long
 * functions for MPI codes than I am for some other settings, but since
 * I have a little more time now to clean things up, let's be good about
 * splitting the code into smaller pieces.  We'll keep a separate struct
 * to manage all the communications.
 *
 * ### Core computation
 *
 * We start with the computational kernel, with a $1/r^2$ interaction function
 * (the same type that we would see in electrostatic or gravitational
 * interactions).  We use unit weights on all the interactions.
 *
 */

void compute_interact(int nother, double* other_pts,
                      int nmine,  double* my_pts,
                      double* results)
{
    for (int i = 0; i < nmine; ++i) {
        double fi = 0.0;
        double xi = my_pts[i];
        for (int j = 0; j < nother; ++j) {
            double dij = xi-other_pts[j];
            if (dij != 0)
                fi += 1.0/(dij*dij);
        }
        results[i] += fi;
    }
}

/**
 * ### The main data structure
 *
 * Now let's set up a struct for all the communication buffers and
 * indexing structures that we need.  The `a2a_init` function sets up
 * the indexing structure and allocates all the space needed;
 * the `a2a_free` tears everything down.
 *
 */
typedef struct a2a_sim_t {

    // -- Number of procs and local rank in MPI communicator
    int nproc;
    int rank;

    // -- Offsets and counts of points owned by each rank
    int* restrict offsets;        // Start offsets for each rank in global idx
    int* restrict counts;         // Counts of data owned at each rank
    int max_count;                // Max counts owned by any rank

    // -- Accumulators and points (local)
    double* restrict results;     // Accumulators for local results
    double* restrict points;      // Local points

    // -- Buffers for pass-it-around protocol
    double* restrict rpoints;     // Remote points for current computation
    double* restrict sendbuf;     // Buffer for sending data
    double* restrict recvbuf;     // Buffer for receiving data
    MPI_Request req[2];           // Send and receive requests

} a2a_sim_t;


// Helper function: check an allocation and clear the space
void* clear_malloc(size_t len)
{
    void* data = malloc(len);
    assert(data != NULL);
    return memset(data, 0, len);
}


// Initialize an a2a_sim_t data structure
void a2a_init(a2a_sim_t* s, int npoints)
{
    // Get size and rank info from communicator
    int nproc;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Set up offsets array
    int* offsets = (int*) malloc((nproc+1) * sizeof(int));
    int* counts = (int*) malloc(nproc * sizeof(int));
    for (int i = 0; i <= nproc; ++i)
        offsets[i] = (i * npoints)/nproc;

    // Compute counts from offsets
    int max_count = 0;
    for (int i = 0; i < nproc; ++i) {
        counts[i] = offsets[i+1]-offsets[i];
        if (counts[i] > max_count)
            max_count = counts[i];
    }

    // Set fields and allocate data
    s->nproc = nproc;
    s->rank  = rank;
    s->offsets = offsets;
    s->counts  = counts;
    s->max_count = max_count;
    s->results = (double*) clear_malloc(counts[rank] * sizeof(double));
    s->points  = (double*) clear_malloc(counts[rank] * sizeof(double));
    s->rpoints = (double*) clear_malloc(max_count * sizeof(double));
    s->sendbuf = (double*) clear_malloc(max_count * sizeof(double));
    s->recvbuf = (double*) clear_malloc(max_count * sizeof(double));
}


// Free resources for an a2a_sim_t
void a2a_free(a2a_sim_t* s)
{
    free(s->recvbuf);
    free(s->sendbuf);
    free(s->rpoints);
    free(s->points);
    free(s->results);
    free(s->counts);
    free(s->offsets);
}


/**
 * ### Scattering and gathering
 *
 * In order to avoid the complexities of parallel random number generation,
 * we are going to create all the initial points at rank 0 and then *scatter*
 * them to the other processors.
 *
 */
void a2a_generate_points(a2a_sim_t* s)
{
    // Initialize the global points array
    int ntotal = s->offsets[s->nproc];
    double* all_points = NULL;
    if (s->rank == 0) {
        all_points = (double*) malloc(ntotal * sizeof(double));
        printf("Initialize array:");
        for (int i = 0; i < ntotal; ++i) {
            all_points[i] = drand48();
            printf(" %g", all_points[i]);
        }
        printf("\n");
    }

    // Distribute points to all ranks (scatter)
    MPI_Scatterv(all_points, s->counts, s->offsets, MPI_DOUBLE,
                 s->points, s->counts[s->rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    // Free the global array
    if (s->rank == 0)
        free(all_points);

    // Print what we got at each rank
    printf("On rank %d (%d-%d): ", s->rank,
           s->offsets[s->rank], s->offsets[s->rank+1]);
    for (int i = 0; i < s->counts[s->rank]; ++i)
        printf(" %g", s->points[i]);
    printf("\n");
}

/**
 * Conversely, at the end we will *gather* all the results to rank 0
 * for output.
 *
 */
void a2a_print_results(a2a_sim_t* s)
{
    // Allocate space to receive results
    int ntotal = s->offsets[s->nproc];
    double* all_results = NULL;
    if (s->rank == 0)
        all_results = malloc(ntotal * sizeof(double));

    // Gather data from all ranks to print out
    MPI_Gatherv(s->results, s->counts[s->rank], MPI_DOUBLE,
                all_results, s->counts, s->offsets, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Print results at rank 0 and free buffer
    if (s->rank == 0) {
        printf("Results: ");
        for (int i = 0; i < ntotal; ++i)
            printf(" %g", all_results[i]);
        printf("\n");
        free(all_results);
    }
}

/**
 * ### Passing data
 *
 * In the main loop, we'll send from each processor to the next processor
 * mod the total number of processors.  We do this with a nonblocking
 * send/receive operation.  Our communication thus has two phases:
 *
 * 1. We copy data into the send buffer and start the send/receive pair.
 * 2. We wait on message completions and copy data out of the recv buffer.
 *
 * In between these two phases, we can do computations that do not involve
 * either of the buffers being communicated.
 *
 */
void a2a_start_sendrecv(a2a_sim_t* s, int phase)
{
    // Previous and next processors in a ring
    int rprev = (s->rank + s->nproc - 1) % s->nproc;
    int rnext = (s->rank + 1) % s->nproc;

    // Copy current remote points into sendbuf space
    memcpy(s->sendbuf, s->rpoints, s->max_count * sizeof(double));

    // Start nonblocking send/recieve
    MPI_Isend(s->sendbuf, s->max_count, MPI_DOUBLE, rnext, phase,
              MPI_COMM_WORLD, s->req+0);
    MPI_Irecv(s->recvbuf, s->max_count, MPI_DOUBLE, rprev, phase,
              MPI_COMM_WORLD, s->req+1);
}


void a2a_end_sendrecv(a2a_sim_t* s)
{
    // Wait for requests to complete
    MPI_Waitall(2, s->req, MPI_STATUSES_IGNORE);

    // Copy data from recvbuf into rpoints for next phase of computation
    memcpy(s->rpoints, s->recvbuf, s->max_count * sizeof(double));
}


/**
 * ### The main computation
 *
 * Having wrapped all the MPI calls in previous helper functions,
 * the main routine to generate points, do all the interactions,
 * and gather and print results can now be fairly terse.
 *
 */
void compute_interacts()
{
    int npoints = 5;
    a2a_sim_t sim;
    a2a_init(&sim, npoints);
    a2a_generate_points(&sim);

    a2a_start_sendrecv(&sim, 0);
    for (int phase = 0; phase < sim.nproc; ++phase) {
        int rdata = ((sim.rank-phase) + sim.nproc) % sim.nproc;
        compute_interact(sim.counts[rdata], sim.rpoints,
                         sim.counts[sim.rank], sim.points, sim.results);
        a2a_end_sendrecv(&sim);
        if (phase < sim.nproc-1)
            a2a_start_sendrecv(&sim, phase+1);
    }

    a2a_print_results(&sim);
    a2a_free(&sim);
}


/**
 * ## The main routine
 */
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
    compute_interacts();

    MPI_Finalize();
    return 0;
}
