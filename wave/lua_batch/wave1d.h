#ifndef WAVE1D_H
#define WAVE1D_H

//ldoc on
/**
 * ## Wave timestepper interface
 * 
 */
typedef struct {
    int n;                       // Number of mesh cells, including boundary
    float c;                     // Wave speed
    float dt;                    // Time step
    float C2;                    // Squared CFL constant (c*dt/dx)^2
    int tidx;                    // Index of most recent time step
    int npart;                   // Number of partitions
    int* restrict offsets;       // Offsets of partitions
    float* restrict us;          // Global data
    float* restrict subdomains;  // Subdomains for fast stepping
} wave1d_t;

/**
 * The `wave1d_sim_t` struct contains the simulation state.  This is a dynamic
 * object that must be allocated by `wave1d_alloc` and freed by `wave1d_free`.
 * The `wave1d_alloc` function also takes some basic parameters:
 * the number of mesh cells `n`, the wave speed `c`, and the time step `dt`.
 *
 */
wave1d_t* wave1d_alloc(int n, float c, float dt);
void wave1d_free(wave1d_t* sim);

/**
 * The time `wave1d_steps` function advances the state of the simulation by 
 * `nsteps` and returns the identifier of the last step taken.  We assume
 * initial conditions are given at steps -1 and 0, respectively.  The
 * `wave1d_steps` function returns the index of the current step.
 * The `wave1d_batch` function returns a natural batch size for advancing
 * steps.
 *
 */
int wave1d_steps(wave1d_t* sim, int nsteps);
int wave1d_batch();

/**
 * We need a way to access the first couple frames (at the very least);
 * we do this with the `wave1d_frame` routine.
 * Finally, `wave1d_dump` writes the simulation state to a file.
 *
 */
float* wave1d_frame(wave1d_t* sim, int step);
void wave1d_dump(wave1d_t* sim, const char* fname);

//ldoc off
#endif // WAVE1D_H
