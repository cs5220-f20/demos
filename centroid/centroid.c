#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

void fill_array(double* p, int n)
{
    for (int i = 0; i < n; ++i)
        p[i] = (double) rand() / RAND_MAX; 
}

void centroid1(double* result, double* xy, int n)
{
    double x = 0.0;
    double y = 0.0;
    for (int i = 0; i < n; ++i) {
        x += xy[2*i+0];
        y += xy[2*i+1];
    }
    result[0] = x/n;
    result[1] = y/n;
}

void centroid2(double* result, double* xy, int n)
{
    double x = 0.0;
    double y = 0.0;
    for (int i = 0; i < n; ++i)
        x += xy[2*i+0];
    for (int i = 0; i < n; ++i)
        y += xy[2*i+1];
    result[0] = x/n;
    result[1] = y/n;
}

void centroid3(double* result, double* xs, double* ys, int n)
{
    double x = 0.0;
    double y = 0.0;
    for (int i = 0; i < n; ++i)
        x += xs[i];
    for (int i = 0; i < n; ++i)
        y += ys[i];
    result[0] = x/n;
    result[1] = y/n;
}

int main()
{
    int n = 1000000;
    double* xy = (double*) calloc(2 * n, sizeof(double));
    double result[2];
    fill_array(xy, 2*n);

    // Timing first option
    double t = omp_get_wtime();
    for (int k = 0; k < 100; ++k)
        centroid1(result, xy, n);
    printf("%g (%g, %g)\n", (omp_get_wtime()-t)/100, result[0], result[1]);

    // Timing second option
    t = omp_get_wtime();
    for (int k = 0; k < 100; ++k)
        centroid2(result, xy, n);
    printf("%g (%g, %g)\n", (omp_get_wtime()-t)/100, result[0], result[1]);

    // Timing third option
    t = omp_get_wtime();
    for (int k = 0; k < 100; ++k)
        centroid3(result, xy, xy+n, n);
    printf("%g (%g, %g)\n", (omp_get_wtime()-t)/100, result[0], result[1]);

    free(xy);
}
