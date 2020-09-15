# Centroid

Three versions of a centroid computation among a million data points.  Try the following:

1.  Run the code as it is (use `make` to build the code).  What do you observe?
2.  What happens if you compile with `-O3` (changing the `CFLAGS` in the Makefile) and
    eliminate the printout of the centroid result?
3.  Try to get the information from GCC about vectorization.  What do you observe?
4.  Try some of the tricks mentioned in the code optimization lecture (aligning the memory,
    explicitly re-associating to reveal independent work, etc).  Can you make things much faster?