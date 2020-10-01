/*
 * Take one time step of a standard finite difference solver for the wave
 * equation
 * $$
 *   u_tt = c u_xx
 * $$
 * on a 1D domain with provided Dirichlet data.
 */
void time_step(int n,             // Total number of grid points
               int b,             // Number of boundary cells at each end
               const double* u0,  // Data two time steps ago
               const double* u1,  // Data one time step ago
               double* u2,        // Data to be computed
               double c,          // Wave speed
               double dx,         // Spatial step size
               double dt)         // Time step size
{
    double C = c*(dt/dx);
    double C2 = C*C;
    for (int i = b; i < n-b; ++i)
        u2[i] = 2*u1[i]-u0[i] + C2*(u1[i-1]-2*u1[i]+u1[i+1]);
}
