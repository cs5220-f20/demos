//ldoc on
/**
 * % OpenMP demos
 *
 * This file has a number of small demos of OpenMP constructs.
 * It must be compiled with the OpenMP flags turned on (see
 * the various `Makefile.in` files to see how to do this with
 * different compilers).
 *
 * Any code that uses OpenMP is going to include the `omp.h` header,
 * where the various OpenMP library routines are declared.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/**
 * ## Hello, world
 * 
 * This routine uses `omp_get_thread_num` and `omp_get_num_threads` to tell
 * us the index of the current thread and the total number of threads in
 * the region
 *
 */
void hello()
{
    #pragma omp parallel
    printf("Hello world from %d / %d\n", 
           omp_get_thread_num(), omp_get_num_threads());
}

/**
 * ## Shared and private
 * 
 * The `hello` routine doesn't do anything interesting with memory.
 * Here are two of the simplest possible examples that actually do something,
 * even if it is just filling an array.  We've included a couple of additional
 * tweaks here that aren't in the slides:
 *
 * - We use `omp_get_max_threads()` to determine the size of the
 *   array that the threads write to.  This tells us the maximum number
 *   of threads allowed.
 * - We use `omp_get_num_threads()` to tell how many threads we actually
 *   got.  Note that `omp_get_num_threads()` tells the number of threads
 *   *in the current parallel section*.  So it has to be called from within
 *   the parallel section, not outside.  We use the `lastprivate` variable
 *   declaration to make sure that the `nthreads` variable gets copied out
 *   to the outside world when the parallel region ends.
 *
 * We also use a variable-length array for the output.  This is technically
 * optional in C11 (it was required in C99), but all the common compilers
 * support it, including CLang, Intel, and GCC.
 *
 */

void print_array(const char* tag, int* array, int n)
{
    printf("%s:", tag);
    for (int i = 0; i < n; ++i)
        printf(" %d", array[i]);
    printf("\n");
}

void fill_array1()
{
    int max_threads = omp_get_max_threads();
    int i;
    int s[max_threads];
    memset(s, 0, max_threads * sizeof(int));
    #pragma omp parallel shared(s) private(i)
    {
        i = omp_get_thread_num();
	s[i] = i;
    }
    print_array("fill_array1:", s, max_threads);
}

void fill_array2()
{
    int max_threads = omp_get_max_threads();
    int s[max_threads];
    memset(s, 0, max_threads * sizeof(int));
    #pragma omp parallel
    {
        int i = omp_get_thread_num();
	s[i] = i;
    }
    print_array("fill_array1:", s, max_threads);
}

/**
 * ## Monte Carlo and reductions
 * 
 * This is a Monte Carlo experiment to compute $\pi/4$.  We are going to
 * run this on several different threads, so we need to use different
 * pseudo-random number streams for each thread.  The `rand48` family of
 * random number generators is not great, but the better alternatives
 * mostly involve bringing in external libraries, so we'll leave be.
 *
 */
double mc_pi(int ntrials, uint64_t seed)
{
    double sumX = 0.0;
    for (int i = 0; i < ntrials; ++i) {
        double x = erand48((unsigned short*) &seed);
        double y = erand48((unsigned short*) &seed);
        if (x*x + y*y < 1)
            sumX += 1.0;
    }
    return sumX;
}

void omp_mc_pi(int ntrials)
{
    int max_threads = omp_get_max_threads();
    uint64_t seeds[max_threads];
    for (int i = 0; i < max_threads; ++i)
        seeds[i] = lrand48();

    double sum = 0.0;
    int all_trials = 0;
    #pragma omp parallel reduction(+:sum) reduction(+:all_trials)
    {
        sum = mc_pi(ntrials, seeds[omp_get_thread_num()]);
        all_trials = ntrials;
    }

    printf("Estimate %g (vs %g) from %d trials\n",
           4*sum/all_trials, M_PI, all_trials);
}

/**
 * ## Linked list and critical sections
 *
 * Critical sections are good for protecting data structures from
 * inconsistencies during updates.  As an example, let's consider a
 * linked list data structure that is concurrently updated by several threads.
 * Note that if we have several linked lists, we still end up with just one
 * critical section!  So finer-grain control might make sense in that case.
 *
 */
typedef struct link_t {
    struct link_t* next;
    int data;
} link_t;

void list_push(link_t** l, int data)
{
    link_t* link = (link_t*) malloc(sizeof(link_t));
    link->data = data;
    #pragma omp critical(list_cs)
    {
        link->next = *l;
        *l = link;
    }
}

void print_list(const char* s, link_t* l)
{
    printf("%s:", s);
    while (l != NULL) {
        printf(" %d", l->data);
        l = l->next;
    }
    printf("\n");
}

void free_list(link_t* l)
{
    while (l != NULL) {
        link_t* lnext = l->next;
        free(l);
        l = l->next;
    }
}

void list_test()
{
    link_t* l = NULL;
    #pragma omp parallel for
    for (int i = 0; i < 10; ++i)
        list_push(&l, i);
    print_list("list_test", l);
    free_list(l);
}

void list_push2(link_t** l, int data)
{
    link_t* link = (link_t*) malloc(sizeof(link_t));
    link->data = data;
    #pragma omp atomic capture
    {
        link->next = *l;
        *l = link;
    }
}

/**
 * ## Linked list and atomic capture
 *
 * The linked list example above uses a critical section to protect
 * the list push operation.  The version below uses an atomic capture
 * operation (i.e. an atomic operation that reads out the old value
 * from some memory and replaces it with a new value).  This atomic
 * capture approach is finer-grained than the critical section; for
 * example, we would not have any contention when using this approach
 * to update two distinct lists concurrently.
 *
 */

void list_test2()
{
    link_t* l = NULL;
    #pragma omp parallel for
    for (int i = 0; i < 10; ++i)
        list_push(&l, i);
    print_list("list_test", l);
    free_list(l);
}

/**
 * ## Atomic updates
 *
 * The mo common use of atomic operations is to update small shared data
 * items - counters, most often.
 *
 */
void count_evens()
{
    int counts[2] = {0, 0};
    #pragma omp parallel for
    for (int i = 0; i < 100; ++i) {
        int which_counter = (i % 2);
        #pragma omp atomic
        counts[which_counter]++;
    }
    printf("For 0 to 99: %d evens, %d odds\n", counts[0], counts[1]);
}

/**
 * ## Barrier-based synchronization
 *
 * Barrier-based synchronization is particularly useful when we have
 * computations organized into distinct time steps or phases, with
 * each step depending on data written in the previous step.  We've seen
 * a few such computations, including the 1D wave equation.  For this
 * example, without attempting to make anything go fast, let's try
 * the Game of Life. 
 *
 * We start with some (non-parallel) setup.  We'll use a simple glider
 * pattern to test things out - this pattern translates by one cell
 * horizontally and vertically once every four generations.
 *
 */
#define NLIFE 100

typedef struct life_board_t {
    uint8_t cells[2][NLIFE][NLIFE];
} life_board_t;

const char* glider =
    " *\n"
    "  *\n"
    "***\n";

void init_life(life_board_t* board, const char* s)
{
    memset(board->cells, 0, 2 * NLIFE * NLIFE);
    for (int i = 1; i < NLIFE-1 && *s != '\0'; ++i) {
        int live = 0;
        for (int j = 1; *s && *s != '\n'; ++j, ++s) {
            if (i < NLIFE-1) {
                if (*s == ' ')
                    board->cells[0][i][j] = 0;
                else if (*s == '*') {
                    board->cells[0][i][j] = 1;
                    ++live;
                }
            }
        }
        if (*s == '\n')
            ++s;
    }
}

void print_life(life_board_t* board, int gen, int nrow, int ncol)
{
    for (int i = 1; i <= nrow; ++i) {
        for (int j = 1; j <= ncol; ++j)
            printf(board->cells[gen][i][j] ? "*" : " ");
        printf("\n");
    }
}

/**
 * The actual `run_life` routine updates the board generation by
 * generation.  At each step, we are only reading from the current
 * generation, and have independent writes into the board for the next
 * generation.  Note that the outer loop (over generations) is *not*
 * parallel: each thread runs all the iterations of this outer loop.
 * But we do decorate the inner loop next to be parallelized across the
 * threads, so each thread in the team updates some part of the board.
 * There is an implicit barrier at the end of a parallel for loop
 * unless we include the `nowait` clause, so each loop switches generations
 * in lockstep.
 *
 */
int run_life(life_board_t* board, int ngens)
{
    #pragma omp parallel
    for (int g = 0; g < ngens; ++g) {
        int current_gen = g%2;
        int next_gen = (g+1)%2;

        #pragma omp for collapse(2)
        for (int i = 1; i < NLIFE-1; ++i) {
            for (int j = 1; j < NLIFE-1; ++j) {

                // Count live neighbors
                int count = -board->cells[current_gen][i][j];;
                for (int ii = -1; ii <= 1; ++ii)
                    for (int jj = -1; jj <= 1; ++jj)
                        count += board->cells[current_gen][i+ii][j+jj];

                // Update rule
                if (board->cells[current_gen][i][j] &&
                    (count == 2 || count == 3))
                    board->cells[next_gen][i][j] = 1; // Still alive
                else if (!board->cells[current_gen][i][j] && count == 3)
                    board->cells[next_gen][i][j] = 1; // Birth
                else
                    board->cells[next_gen][i][j] = 0; // Dead
            }
        } // Implicit barrier at end of parallel for
    }
    return ngens%2;
}

void glider_test()
{
    life_board_t board;
    init_life(&board, glider);
    printf("... Gen 0... :\n");
    print_life(&board, 0, 6, 8);
    run_life(&board, 12);
    printf("... Gen 12... :\n");
    print_life(&board, 0, 6, 8);
}

/**
 * ## Task-based list traversal
 *
 * This list traversal code goes through a linked list and creates a task
 * for processing each.
 *
 */
void list_traversal_test()
{
    link_t* head = NULL;
    for (int i = 0; i < 10; ++i)
        list_push2(&head, i);

    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            printf("Creating tasks from %d\n", omp_get_thread_num());
            for (link_t* link = head; link; link = link->next) {
                #pragma omp task firstprivate(link)
                printf("Process %d on thread %d\n",
                       link->data, omp_get_thread_num());
            }
            printf("Done creating tasks\n");
         }
     }
}

// ldoc off
int main()
{
    printf("\n--- hello ---\n");
    hello();

    printf("\n--- fill array --\n");
    fill_array1();
    fill_array2();

    printf("\n--- Monte Carlo pi --\n");
    omp_mc_pi(10000);

    printf("\n--- List and critical section --\n");
    list_test();

    printf("\n--- List and atomic --\n");
    list_test2();

    printf("\n--- Atomic counter updates --\n");
    count_evens();

    printf("\n--- Game of Life --\n");
    glider_test();

    printf("\n--- List traversal --\n");
    list_traversal_test();
}
