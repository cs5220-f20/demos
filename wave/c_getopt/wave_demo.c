#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


//ldoc on
/**
 * % The 1D wave equation (C edition)
 *
 * This is another serial implementation of the wave stepper code.
 * The purpose is to give some tips about how to organize a time-stepper
 * code like this (e.g. using a limited history buffer for time steps),
 * and also to show off how the argument processing works.
 *
 * This is also an opportunity for me to show of my
 * [`ldoc`](https://github.com/dbindel/ldoc) documentation generator.
 *
 * ## Argument processing
 *
 * The [`getopt`](https://en.wikipedia.org/wiki/Getopt)
 * command is a parser for Unix-style command line flags.
 * It's a bit clunky and awkward, but still worth knowing about.
 * Even a code as simple as the 1D wave equation solver has several
 * possible parameters: the number of grid points, the wave speed,
 * the time step, the number of time steps, and the name of the output
 * file, in this case.
 *
 */
int parse_args(int argc, char** argv,
               int* n, double* c, double* dt, int* nsteps, char** ofname)
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

/**
 * ## The time stepper
 *
 * The `time_step` command takes one time step of a standard finite
 * difference solver for the wave equation
 * $$
 *   u_{tt} = c^2 u_{xx}
 * $$
 * on a 1D domain with provided Dirichlet data.  We approximate
 * both second derivatives by finite differences, i.e.
 * $$
 *   u_{tt}(x,t) \approx
 *     \frac{u(x,t-\Delta t) - 2 u(x,t) + u(x,t+\Delta t)}{(\Delta t)^2}
 * $$
 * and similarly for the second derivative in space.  Plugging this
 * into the wave equation and rearranging terms gives us the update
 * formula
 * $$
 *   u(x, t + \Delta t) =
 *     2 u(x,t) - u(x,t-\Delta t) +
 *     \left( c \frac{\Delta t}{\Delta x} \right)^2
 *     \left( u(x+\Delta x, t) - 2 u(x,t) + u(x-\Delta x, t) \right)
 * $$
 * We let `u0[i]` and `u1[i]` denote the values of $u(x_i,t-\Delta t)$ and
 * $u(x_i, t)$ for a mesh of points $x_i$ with spacing $\Delta x$,
 * and write to the output `u2` which contains $u(x, t + \Delta t)$.
 *
 * The wave equation requires boundary conditions as well as initial
 * data, and this time stepper does not explicitly handle boundary conditions.
 * Instead, the boundary conditions need to be handled by a separate routine
 * that updates a layer of `b` "ghost cells" at each end of the computational
 * mesh.  This can be done by filling in fixed values (for Dirichlet boundary
 * conditions), copying endpoint values (for homogeneous Neumann conditions),
 * or transferring data from a neighboring processor (for the case where the
 * time stepper is being used inside of another code).
 *
 */
void time_step(int n,             // Total number of grid points
               int b,             // Number of boundary cells at each end
               const double* restrict u0,  // Data two time steps ago
               const double* restrict u1,  // Data one time step ago
               double* restrict u2,        // Data to be computed
               double c,          // Wave speed
               double dx,         // Spatial step size
               double dt)         // Time step size
{
    double C = c*(dt/dx);
    double C2 = C*C;
    for (int i = b; i < n-b; ++i)
        u2[i] = 2*u1[i]-u0[i] + C2*(u1[i-1]-2*u1[i]+u1[i+1]);
}

/**
 * Of course, we generally don't want to take exactly one time step:
 * we want to take many time steps.  We only need to keep three
 * consecutive time steps.  We indicate the most recent completed step
 * with `i0`; the step before is at `i0-1` (mod 3) and the step to be
 * computed is at `i0+1` (mod 3).  We'll return the index (mod 3)
 * of the last time step computed.
 *
 */
int time_steps(int n, int b, double* us, int i0, int nsteps,
               double c, double dx, double dt)
{
    for (int j = 0; j < nsteps; ++j) {
        double* u0 = us + ((i0+2+j)%3)*n;
        double* u1 = us + ((i0+0+j)%3)*n;
        double* u2 = us + ((i0+1+j)%3)*n;
        time_step(n, 1, u0, u1, u2, c, dx, dt);
    }
    return (i0+nsteps)%3;
}

/**
 * ## Initial conditions
 *
 * If we want a general purpose code, it often makes sense to separate
 * out the initial conditions from everything else.  Certainly you
 * want to make it easy to swap out one set of initial conditions for
 * another.
 *
 */
void initial_conditions(int n, double* u0, double* u1,
                        double c, double dx, double dt)
{
    for (int i = 1; i < n-1; ++i) {
        double x = i*dx;
        u0[i] = exp(-25*(x-0.5)*(x-0.5));
        u1[i] = exp(-25*(x-0.5-c*dt)*(x-0.5-c*dt));
    }
}

/**
 * ## Printing the state
 *
 * We will use a simple text file format for the simulation state.
 * That makes it easier to read into a viewer like the Python script
 * that we wrote for the Python version of this code.  A binary format
 * (or a compressed binary format) might be worth considering if we are
 * planning to output more data.  One might also consider downsampling
 * the output (e.g. printing only every mth mesh point) if the purpose
 * of the output is visualization rather than detailed diagnosis
 * or restarting of a simulation.
 *
 */
void print_mesh(FILE* fp, int n, double* u)
{
    for (int i = 0; i < n; ++i)
        fprintf(fp, "%g\n", u[i]);
}


/**
 * ## The `main` event
 *
 */
int main(int argc, char** argv)
{
    // Set defaults and parse arguments
    int n = 1000;
    double c = 1.0;
    double dt = 5e-4;
    int nsteps = 3600;
    char* fname = NULL;
    int flag = parse_args(argc, argv, &n, &c, &dt, &nsteps, &fname);
    if (flag != 0)
        return flag;

    // Compute space step and check CFL
    double dx = 1.0/(n-1);
    double C = c*dt/dx;
    if (C >= 1.0) {
        printf("CFL constant is %g (should be < 1 for stability)\n", C);
        return -1;
    }

    // Setting up the storage space for the time steps
    double* us = (double*) malloc(3 * n * sizeof(double));
    memset(us, 0, 3 * n * sizeof(double));

    // Initialize the problem and run the time stepper
    initial_conditions(n, us+0*n, us+1*n, c, dx, dt);
    int last_idx = time_steps(n, 1, us, 1, nsteps, c, dx, dt);

    // Write the final output
    if (fname) {
        FILE* fp = fopen(fname, "w");
        print_mesh(fp, n, us+n*last_idx);
        fclose(fp);
    }

    // Clean up and return
    free(us);
    return 0;
}
