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

void centroid1a(double* restrict result, 
                const double* restrict xy, int n)
{
    double x = 0.0;
    double y = 0.0;
    #pragma omp simd reduction(+:x) reduction(+:y)
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

void centroid3a(double* restrict result,
                const double* restrict xs,
                const double* restrict ys,
                int n)
{
    double x = 0.0;
    double y = 0.0;

    #pragma omp simd reduction(+:x)
    for (int i = 0; i < n; ++i)
        x += xs[i];

    #pragma omp simd reduction(+:x)
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
    printf("1:  %g (%g, %g)\n", (omp_get_wtime()-t)/100, result[0], result[1]);

   // Timing first option
    t = omp_get_wtime();
    for (int k = 0; k < 100; ++k)
        centroid1a(result, xy, n);
    printf("1a: %g (%g, %g)\n", (omp_get_wtime()-t)/100, result[0], result[1]);

    // Timing second option
    t = omp_get_wtime();
    for (int k = 0; k < 100; ++k)
        centroid2(result, xy, n);
    printf("2:  %g (%g, %g)\n", (omp_get_wtime()-t)/100, result[0], result[1]);

    // Timing third option
    t = omp_get_wtime();
    for (int k = 0; k < 100; ++k)
        centroid3(result, xy, xy+n, n);
    printf("3:  %g (%g, %g)\n", (omp_get_wtime()-t)/100, result[0], result[1]);

    // Timing third option
    t = omp_get_wtime();
    for (int k = 0; k < 100; ++k)
        centroid3a(result, xy, xy+n, n);
    printf("3a: %g (%g, %g)\n", (omp_get_wtime()-t)/100, result[0], result[1]);

    // Compare "speed of light"
    double mem_bw = 12.4e9;  // GB/s
    double mem_read = 2*n*sizeof(double);
    printf("Time to read data: %g\n", mem_read/mem_bw);

    free(xy);
}
