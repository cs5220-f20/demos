#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "waves.h"

// Demonstrate the use of getopt for argument parsing
int parse_args(int argc, char** argv, int* n, double* c, double* dt, int* nsteps, char** ofname)
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
        "  o: Output file name (default '%s')\n";

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


void print_mesh(FILE* fp, int n, double* u)
{
    for (int i = 0; i < n; ++i)
        fprintf(fp, "%g\n", u[i]);
}


int main(int argc, char** argv)
{
    // Set defaults and parse arguments
    int n = 1000;
    double c = 1.0;
    double dt = 5e-4;
    int nsteps = 3600;
    char* fname = "steps.txt";
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
    double* u0 = us + 0*n;
    double* u1 = us + 1*n;
    double* u2 = us + 2*n;
    
    // Set up initial conditions
    for (int i = 1; i < n-1; ++i) {
        double x = i*dx;
        u0[i] = exp(-25*(x-0.5)*(x-0.5));
        u1[i] = exp(-25*(x-0.5-c*dt)*(x-0.5-c*dt));
    }
    
    // Run the time stepper
    for (int j = 0; j < nsteps; ++j) {
        u0 = us + ((0+j)%3)*n;
        u1 = us + ((1+j)%3)*n;
        u2 = us + ((2+j)%3)*n;
        time_step(n, 1, u0, u1, u2, c, dx, dt);
    }
    
    // Write the final output
    FILE* fp = fopen(fname, "w");
    print_mesh(fp, n, us);
    fclose(fp);
    
    free(us);
    return 0;
}