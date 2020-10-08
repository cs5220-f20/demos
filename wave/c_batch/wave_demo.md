% The 1D wave equation (C edition)

This is another serial implementation of the wave stepper code.
Relative to the earlier C version, this one is tuned a bit to
try to get better (serial) performance.

This is also an opportunity for me to show of my
[`ldoc`](https://github.com/dbindel/ldoc) documentation generator.

## Argument processing

The [`getopt`](https://en.wikipedia.org/wiki/Getopt)
command is a parser for Unix-style command line flags.
It's a bit clunky and awkward, but still worth knowing about.
Even a code as simple as the 1D wave equation solver has several
possible parameters: the number of grid points, the wave speed,
the time step, the number of time steps, and the name of the output
file, in this case.

    int parse_args(int argc, char** argv,
                   int* n, float* c, float* dt, int* nsteps, char** ofname)
    {
        int opt;
        int err_flag = 0;
        const char* help_string =
            "Usage: wave_demo [-h] [-n nmesh] [-c speed] [-t dt] [-s nsteps] [-o fname]\n"
            "\n"
            "  h: Print this message\n"
            "  n: Set the mesh size (default %d)\n"
            "  c: Set the speed of sound (default %g)\n"
            "  t: Set the time step (default %g)\n"
            "  s: Set the number of steps (default %d)\n"
            "  o: Output file name (default none)\n";
    
        while ((opt = getopt(argc, argv, "hn:c:t:s:o:")) != -1) {
            switch (opt) {
                case 'h':
                    fprintf(stderr, help_string, *n, *c, *dt, *nsteps, *ofname);
                    err_flag = 1;
                    break;
                case 'n':
                    *n = atoi(optarg);
                    if (*n < 3) {
                        fprintf(stderr, "Error: Need at least three mesh points\n");
                        err_flag = -1;
                    }
                    break;
                case 'c':
                    *c = atof(optarg);
                    break;
                case 't':
                    if (*dt <= 0.0) {
                        fprintf(stderr, "Error: Time step must be positive\n");
                        err_flag = 1;
                    }
                    *dt = atof(optarg);
                    break;
                case 's':
                    *nsteps = atoi(optarg);
                    break;
                case 'o':
                    *ofname = strdup(optarg);
                    break;
                case '?':
                    fprintf(stderr, "Unknown option\n");
                    err_flag = -1;
                    break;
            }
        }
        return err_flag;
    }
    
## The time stepper

The `time_step` command takes one time step of a standard finite
difference solver for the wave equation
$$
  u_{tt} = c^2 u_{xx}
$$
on a 1D domain with provided Dirichlet data.  We approximate
both second derivatives by finite differences, i.e.
$$
  u_{tt}(x,t) \approx
    \frac{u(x,t-\Delta t) - 2 u(x,t) + u(x,t+\Delta t)}{(\Delta t)^2}
$$
and similarly for the second derivative in space.  Plugging this
into the wave equation and rearranging terms gives us the update
formula
$$
  u(x, t + \Delta t) =
    2 u(x,t) - u(x,t-\Delta t) +
    \left( c \frac{\Delta t}{\Delta x} \right)^2
    \left( u(x+\Delta x, t) - 2 u(x,t) + u(x-\Delta x, t) \right)
$$
We let `u0[i]` and `u1[i]` denote the values of $u(x_i,t-\Delta t)$ and
$u(x_i, t)$ for a mesh of points $x_i$ with spacing $\Delta x$,
and write to the output `u2` which contains $u(x, t + \Delta t)$.

    void time_step(int n,             // Total grid points to update
                   const float* restrict u0,  // Data two time steps ago
                   const float* restrict u1,  // Data one time step ago
                   float* restrict u2,        // Data to be computed
                   float C2)          // Squared nondim wave speed
    {
        float C1 = 2.0f-2.0f*C2;
        for (int i = 0; i < n; ++i)
            u2[i] = C1*u1[i]-u0[i] + C2*(u1[i-1]+u1[i+1]);
    }
    
Of course, we generally don't want to take exactly one time step:
we want to take many time steps.  We only need to keep three
consecutive time steps.  We indicate the most recent completed step
with `i0`; the step before is at `i0-1` (mod 3) and the step to be
computed is at `i0+1` (mod 3).  We'll return the index (mod 3)
of the last time step computed.

We distinguish between the total number of cells that we are keeping
at each step (`ntotal`) and the number of cells that we are updating
(`nupdate`).  In general, something is wrong if the ntotal < nupdate+2.
    int time_steps(int ntotal,   // Total number of cells (including ghosts)
                   int nupdate,  // Number of cells to update
                   float* us,   // Start of cells to be updated
                   int i0,       // Index of most recent completed step
                   int nsteps,   // Number of steps to advance
                   float C2)    // Squared wave speed
    {
        assert(ntotal >= nupdate+2);
        for (int j = 0; j < nsteps; ++j) {
            float* u0 = us + ((i0+2+j)%3)*ntotal;
            float* u1 = us + ((i0+0+j)%3)*ntotal;
            float* u2 = us + ((i0+1+j)%3)*ntotal;
            time_step(nupdate, u0, u1, u2, C2);
        }
        return (i0+nsteps)%3;
    }
    
## Initial conditions

If we want a general purpose code, it often makes sense to separate
out the initial conditions from everything else.  Certainly you
want to make it easy to swap out one set of initial conditions for
another.

    void initial_conditions(int n, float* u0, float* u1,
                            float c, float dx, float dt)
    {
        for (int i = 1; i < n-1; ++i) {
            float x = i*dx;
            u0[i] = exp(-25*(x-0.5)*(x-0.5));
            u1[i] = exp(-25*(x-0.5-c*dt)*(x-0.5-c*dt));
        }
    }
    
## Printing the state

We will use a simple text file format for the simulation state.
That makes it easier to read into a viewer like the Python script
that we wrote for the Python version of this code.  A binary format
(or a compressed binary format) might be worth considering if we are
planning to output more data.  One might also consider downsampling
the output (e.g. printing only every mth mesh point) if the purpose
of the output is visualization rather than detailed diagnosis
or restarting of a simulation.

    void print_mesh(FILE* fp, int n, float* u)
    {
        for (int i = 0; i < n; ++i)
            fprintf(fp, "%g\n", u[i]);
    }
    
    
## Subdomain partitioning

We break the interior of the domain into disjoint subdomains, each
with index sets `own_start <= i < own_end`.  To advance these owned
mesh points by B steps, we need to look at the range
`sub_start <= i < sub_end` where
`sub_start = min(own_start-B, 0)` and
`sub_end = max(own_end+B, n)`.

We set the default batch size to 40, and the subdomain interior sizes
to 560.

    #define BATCH 40
    #define NS_INNER 1200
    #define NS_TOTAL 1280
    
    inline
    int sub_start(int own_start)
    {
        return (own_start < BATCH) ? 0 : own_start-BATCH;
    }
    
    inline
    int sub_end(int own_end, int n)
    {
        return (own_end+BATCH > n) ? n : own_end+BATCH;
    }
    
The `sub_copyin` routine moves the range from `sub_start` to `sub_end`
into a local copy.  The `sub_copyout` routine moves the range corresponding
to `own_start` to `own_end` (starting at offset `own_start-sub_start`)
from the local array back into the global array.

    void sub_copyin(float* restrict ulocal,
                    float* restrict uglobal,
                    int own_start, int own_end, int n)
    {
        int dstart = sub_start(own_start);
        int dend = sub_end(own_end, n);
        for (int j = 0; j < 3; ++j)
            memcpy(ulocal + j*NS_TOTAL,
                   uglobal + dstart + j*n,
                   (dend - dstart) * sizeof(float));
    }
    
    void sub_copyout(float* restrict ulocal,
                     float* restrict uglobal,
                     int own_start, int own_end, int n)
    {
        int dstart = sub_start(own_start);
        int dend = sub_end(own_end, n);
        for (int j = 0; j < 3; ++j)
            memcpy(uglobal + own_start + j*n,
                   ulocal + own_start - dstart + j*NS_TOTAL,
                   (own_end - own_start) * sizeof(float));
    }
    
The `sub_steps` routine copies data into a local buffer, then advances
that buffer by the given number of steps.
    int sub_steps(int own_start,   // Start of owned cell range
                  int own_end,     // End of owned cell range
                  float* uglobal, // Global storage
                  int n,           // Total number of mesh cells
                  float* ulocal,  // Start of subdomain storage
                  int i0,          // Index of most recent completed step
                  int nsteps,      // Number of steps to advance
                  float C2)       // Squared wave speed
    {
        int dstart = sub_start(own_start);
        int dend = sub_end(own_end, n);
        sub_copyin(ulocal, uglobal, own_start, own_end, n);
        return time_steps(NS_TOTAL, dend-dstart-2, ulocal+1, i0%3, nsteps, C2);
    }
    
The partitioner cuts the domain into pieces of size at most `NS_INNER`;
partition `j` is indices `offsets[j] <= i < offsets[j+1]`.  We return
the total number of partitions in the output argument `npart`; the
offsets array has length `npart+1`.
    int* alloc_partition(int n, int* npart)
    {
        int np = (n-2 + NS_INNER-1)/NS_INNER;
        int* offsets = (int*) malloc((np+1) * sizeof(int));
        for (int i = 0; i <= np; ++i) {
            long r = i*(n-2);
            offsets[i] = 1 + (int) (r/np);
        }
        *npart = np;
        return offsets;
    }
    
## The `main` event

    int main(int argc, char** argv)
    {
        // Set defaults and parse arguments
        int n = 1000;
        float c = 1.0;
        float dt = 5e-4;
        int nsteps = 3600;
        char* fname = NULL;
        int flag = parse_args(argc, argv, &n, &c, &dt, &nsteps, &fname);
        if (flag != 0)
            return flag;
    
        // Compute space step and check CFL
        float dx = 1.0/(n-1);
        float C = c*dt/dx;
        float C2 = C*C;
        if (C >= 1.0) {
            printf("CFL constant is %g (should be < 1 for stability)\n", C);
            return -1;
        }
    
        // Setting up the storage space and initial conditions
        float* us = (float*) malloc(3 * n * sizeof(float));
        memset(us, 0, 3 * n * sizeof(float));
        initial_conditions(n, us+0*n, us+1*n, c, dx, dt);
    
        // Partition the domain
        int npart;
        int* offsets = alloc_partition(n, &npart);
    
        // Set up storage for subdomains
        float* subdomains = (float*) malloc(6 * NS_TOTAL * sizeof(float));
        memset(subdomains, 0, 6 * NS_TOTAL * sizeof(float));
    
        // Run the time stepper in batches of at most B steps
        double t0 = omp_get_wtime();
        int last_idx = (nsteps+1)%3;
        for (int j = 0; j < nsteps; j += BATCH) {
            int bsteps = (j+BATCH > nsteps) ? nsteps-j : BATCH;
            if (npart == 1)
                time_steps(n, n-2, us+1, (j+1)%3, bsteps, C2);
            else {
                for (int i = 0; i <= npart; ++i) {
                    if (i < npart)
                        sub_steps(offsets[i], offsets[i+1], us, n,
                                  subdomains + (i%2) * 3 * NS_TOTAL,
                                  (j+1)%3, bsteps, C2);
                    if (i > 0)
                        sub_copyout(subdomains + ((i+1)%2) * 3 * NS_TOTAL,
                                    us, offsets[i-1], offsets[i], n);
                }
            }
        }
    
        double tfinal = omp_get_wtime();
        double flops = 5.0 * nsteps * n;
        printf("Elapsed time: %g\n", tfinal-t0);
        printf("Effective Gflop/s: %g\n", flops/(tfinal-t0)/1e9);
    
        // Write the final output
        if (fname) {
            FILE* fp = fopen(fname, "w");
            print_mesh(fp, n, us+n*last_idx);
            fclose(fp);
        }
    
        // Clean up and return
        free(subdomains);
        free(offsets);
        free(us);
        return 0;
    }

