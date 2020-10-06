#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "waves.h"

void print_mesh(FILE* fp, int n, double* u)
{
    for (int i = 0; i < n; ++i)
        fprintf(fp, "%g\n", u[i]);
}


int main(int argc, char** argv)
{
    int n = (argc <= 1) ? 1000 : atoi(argv[1]);
    double c = (argc <= 2) ? 1.0 : atof(argv[2]);
    double dt = (argc <= 3) ? 5e-4 : atof(argv[3]);
    double dx = 1.0/(n-1);
    double C = c*dt/dx;
    printf("CFL constant is %g (should be < 1 for stability)\n", C);

    // Setting up the storage space for the time steps
    double* us = (double*) malloc(3 * n * sizeof(double));
    memset(us, 0, 3 * n * sizeof(double));
    double* u0 = us + 0*n;
    double* u1 = us + 1*n;
    double* u2 = us + 2*n;
    
    // Set up initial conditions
    printf("Initialize\n");
    for (int i = 1; i < n-1; ++i) {
        double x = i*dx;
        u0[i] = exp(-25*(x-0.5)*(x-0.5));
        u1[i] = exp(-25*(x-0.5-c*dt)*(x-0.5-c*dt));
    }
    
    // Run the time stepper
    printf("Run time stepper\n");
    int nsteps = 3600;
    for (int j = 0; j < nsteps; ++j) {
        u0 = us + ((0+j)%3)*n;
        u1 = us + ((1+j)%3)*n;
        u2 = us + ((2+j)%3)*n;
        time_step(n, 1, u0, u1, u2, c, dx, dt);
    }
    
    // Write the final output
    FILE* fp = fopen("steps.txt", "w");
    print_mesh(fp, n, us);
    fclose(fp);
    
    free(us);
    return 0;
}