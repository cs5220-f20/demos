# Conway's Game of Life

Conway's "Game of Life" is a 2D cellular automata example.  It is a
great game, interesting mathematically, and fun to play with.
And performance tuning of an implementation of the Game of Life is
also mostly a matter of play; this is not a kernel that serves as
a building block to any more serious application that I know about.
Nonetheless, there are elements of this performance tuning that do
bear on other examples, and so playing about may teach us something.

## I/O costs

The update rules for the Game are cheap.  Output is expensive.  If we
are going to plot or print the board at every generation, it is
probably not worth optimizing anything about the computation.

We could instead consider plotting every several generations, or
keeping a history of some statistics (e.g. number of live cells)
per generation, or plotting only when some interesting event happens.
But we'd probably need to do a lot of generations to come close to the
cost of one board plot.

## Dense and dilute

We are going to discuss the "dense" case of the game, where it makes
sense to update every cell.  We do this partly because it makes a good
building block for the "dilute" case.  It's always nice to have big
blocks of well-structured work that we can rearrange; it makes it
easier to amortize the overheads associated with doing something clever.

If we seriously wanted to try to accelerate the dilute case (where
most of the board is dead), we would probably carve the board into
blocks, and keep a marker for each block saying whether something was
alive or not.  Then we would only bother updating a block when it was
previously alive or neighbored a previously-alive block.

If we were particularly ambitious and memory was an issue, we might
want to dynamically allocate blocks only when they were needed (and
keep an array of pointers to blocks, with NULL pointers for blocks
that were all dead).

## Storage format

Each cell can only have two states: alive or dead.  So we only need
one bit per cell.  A 32K cache has 2^18 bits in it, so we can easily
store a 512-by-512 board completely in L1 cache.  Of course, we need
storage for at least two boards, one for the current generation and
one for the next generation; we would ususally switch back and forth
between them during the game.  Maybe we would want both of those to
fit in cache.

In the dilute case, I would be inclined to organize everything around
48-by-48 blocks, for reasons discussed below.

## Vectorization

There are all sorts of fun games we can play with bit-fiddling in
order to effectively vectorize the cell updates.  There is a [nice
discussion
here](https://www.cs.uaf.edu/2012/spring/cs641/lecture/02_14_SIMD.html).
It's worth knowing that you can do this, and sometimes these types of
bit-fiddling hacks really are useful for inner loops.  But they make
your code look like it was written by a Martian.  So make sure that
you write a reference version that is easier for humans to read,
first.

## Blocking

My strategy if I actually cared about large-scale Game of Life
simulations would probably be this:

- Break the board into 48-by-48 tiles.
- For each tile, update in steps of 8 generations:
  - Read the tile + 8 layers in each direction (a 64-by-64 section in
    local storage)
  - Use the vectorized operations to update the 64-by-64 chunk eight
    times.  At each step, one more boundary layer will be "wrong";
    after eight generations, we will have to throw away the outer
    eight layers of our result.
  - Write the inner 48-by-48 section to the new board section.

This takes advantage of the vector updating cleanly, and lets us do
lots of work on small local working sets.  It's probably a good idea
even on a single core, because of the advantages of doing lots on a
small working set.  It's certainly useful if we ever want to
parallelize; synchronization and communication overheads are likely to
dominate the actual computation costs unless we do something like this.
