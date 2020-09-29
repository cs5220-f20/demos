# Centroid

Three versions of a centroid computation among a million data points.  Try the following:

1.  Run the code as it is (use `make` to build the code).  What do you observe?
2.  What happens if you compile with `-O3` (changing the `CFLAGS` in the Makefile) and
    eliminate the printout of the centroid result?
3.  Try to get the information from GCC about vectorization.  What do you observe?
4.  Try some of the tricks mentioned in the code optimization lecture (aligning the memory,
    explicitly re-associating to reveal independent work, etc).  Can you make things much faster?

## Updated centroid

Since some time has elapsed since we started this example, let's bring it to a successful
conclusion.  If we explicitly use the `-ffast-math` flag, we tell GCC that it is free to
reassociate addition, and give it much more flexibility to reassociate things.  We still
do better with the version that keeps two parallel sums simultaneously.  We include an
example of the OpenMP SIMD directives to try to coax the compiler to vectorize the loop,
but it does nothing that GCC isn't willing to do for us already (at least on
the version installed on Debian Buster).  We've updated the code to include a
"speed of light" computation based on the estimated per core memory bandwidth
on this machine; the current version essentially reaches that speed of light.