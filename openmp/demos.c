#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//ldoc on
/**
 * # OpenMP demos
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
#include <omp.h>

//ldoc off

void print_array(const char* tag, int* array, int n)
{
    printf("%s:", tag);
    for (int i = 0; i < n; ++i)
        printf(" %d", array[i]);
    printf("\n");
}


//ldoc on
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
 * - We use a variable-length array for the output.  This is technically
 *   optional in C11 (it was required in C99), but all the common compilers
 *   support it, including CLang, Intel, and GCC.
 */

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

/**
 * ## Linked list and atomic capture
 *
 */
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

void list_test2()
{
    link_t* l = NULL;
    #pragma omp parallel for
    for (int i = 0; i < 10; ++i)
        list_push(&l, i);
    print_list("list_test", l);
    free_list(l);
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
}
