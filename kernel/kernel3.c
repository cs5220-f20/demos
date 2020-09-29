#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdalign.h>
#include <math.h>
#include <omp.h>


#define D 100
typedef struct {
    alignas(32) float vec[D];
} fvec_t;


void fill_dvec(float* v, int n)
{
    for (int i = 0; i < n; ++i)
        v[i] = (float) drand48();
}


void fill_fvecs(fvec_t* v, int n)
{
    for (int i = 0; i < n; ++i)
        fill_dvec(v[i].vec, D);
}


// This could maybe be treated with an OpenMP SIMD pragma
float feature_dist2(const fvec_t* restrict x,
                    const fvec_t* restrict y)
{
    float r = 0.0;
    for (int i = 0; i < D; ++i) {
        float di = x->vec[i]-y->vec[i];
        r += di*di;
    }
    return r;
}


// Compute f[i] = sum_j k(xi, xj) * c[j] for some kernel k
void kernel_sums_diag(int n,
                      const fvec_t* restrict x,
                      const float* restrict c,
                      float* restrict f)
{
    for (int i = 0; i < n; ++i) {
        f[i] += c[i];
        for (int j = 0; j < i; ++j) {
            float r2 = feature_dist2(x+i, x+j);
            float kr = exp(-r2/2);
            f[i] += kr * c[j];
            f[j] += kr * c[i];
        }
    }
}


void kernel_sums_offdiag(int m, int n,
                         const fvec_t* restrict x,
                         const fvec_t* restrict y,
                         const float* restrict cx,
                         const float* restrict cy,
                         float* restrict fx,
                         float* restrict fy)
{
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            float r2 = feature_dist2(x+i, y+j);
            float kr = exp(-r2/2);
            fx[i] += kr * cy[j];
            fy[j] += kr * cx[j];
        }
    }
}


void kernel_sums(int n,
                 const fvec_t* restrict x,
                 const float* restrict c,
                 float* restrict f)
{
    #define B 10
    memset(f, 0, n * sizeof(float));
    for (int i = 0; i < n; i += B) {
        int bi = (n-i < B) ? n-i : B;
        kernel_sums_diag(bi, x+i, c+i, f+i);
        for (int j = 0; j < i; j += B) {
            int bj = (n-j < B) ? n-j : B;
            kernel_sums_offdiag(bi, bj, x+i, x+j, c+i, c+j, f+i, f+j);
        }
    }
}


int main()
{
    int n = 5000;

    // Set up storage and fill with random stuff
    fvec_t* x = (fvec_t*) aligned_alloc(32, n * sizeof(fvec_t));
    float* f = (float*) malloc(n * sizeof(float));
    float* c = (float*) malloc(n * sizeof(float));

    fill_fvecs(x, n);
    fill_dvec(c, n);

    float t = omp_get_wtime();
    kernel_sums(n, x, c, f);
    t = omp_get_wtime()-t;
    printf("Time: %g\n", t);

    float s = 0.0;
    for (int i = 0; i < n; ++i)
        s += f[i];
    printf("  (dummy result: %g)\n", s);

    free(c);
    free(f);
    free(x);
}
