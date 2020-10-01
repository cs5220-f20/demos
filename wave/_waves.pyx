cimport numpy as np
np.import_array()

# C external function declaration
cdef extern from "waves.h":
    void time_step(int n, int b, double* u0, double* u1, double* u2, 
                   double c, double dx, double dt)

# Python wrapper function
def time_step_func(b,
                   np.ndarray[double, ndim=1, mode="c"] u0 not None,
                   np.ndarray[double, ndim=1, mode="c"] u1 not None,
                   np.ndarray[double, ndim=1, mode="c"] u2 not None,
                   c, dx, dt):
    time_step(u0.shape[0], b,
              <double*> np.PyArray_DATA(u0),
              <double*> np.PyArray_DATA(u1),
              <double*> np.PyArray_DATA(u2),
              c, dx, dt)
