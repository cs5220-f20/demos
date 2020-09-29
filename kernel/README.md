# Kernel matvec example

In this example, we illustrate multiplication of a large kernel matrix
times a vector.  The kernel matrix has the form K(i,j) = phi(rij)
where rij is the Euclidean distance between points x(i) and x(j) in
some high-dimensional space.  For concreteness, we will choose 5000 points
in a 100-dimensional space.

We consider seven different executables, each incrementally faster than
the last.  By a series of careful tweaks, we arrive at an
implementation that is almost 40 times faster than the original
version (running experiments on my Macbook Air with GCC 10).  At no
point do we drop down to the level of writing assembly code, or even
using OpenMP SIMD pragmas.

1.  `vec-1a.x`: A straightforward implementation, compiled with no
    special flags (8.3013s).
2.  `vec0a.x`: Like the previous version, but uses aligned
    declarations for the feature vectors in order to make it easier
    for the compiler to vectorize the distance computations (8.1550).
3.  `vec0b.x`: Like the previous version, but compiled with `-O3` and
    `-march=native` (3.0582).
4.  `vec0c.x`: Like the previous version, but additionally compiled
    with `-ffast-math` (1.5062)
5.  `vec1c.x`: Like the previous version, but takes advantage of
    symmetry of the kernel matrix (K(i,j) = K(j,i)) and avoids
    re-computing matrix elements (0.6127).
6.  `vec2c.x`: Like the previous version, but does some memory
    blocking to reduce cache misses (0.4493s).
7.  `vec3c.x`: Like the previous version, but uses single precision
    rather than double precision (0.2467s).

Note that just the cost of computing the squared distances between each
pair of points costs `d*n*(n+1)` flops (about equal numbers of adds and
multiplies).  The most we could possibly hope for is given by the number
of flops over the peak flop rate, which is just below 0.1 (in double
precision) on this machine.  Hence, we might hope to get another five-fold
speedup or so if we continued to fiddle with the code, and managed to
balance everything perfectly.  Of course, that ignores the time to
actually evaluate the kernel function!
