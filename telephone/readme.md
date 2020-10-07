# Telephone with MPI
-------------------

In this HW, we will cover some Message Passing Interface (MPI) basics. If you have any questions, don't hesitate to post on Piazza or come to OH! 

**Submission:** A tarball/git-link of the provided directory

# The Game of Telephone, with Children

In the children's game of telephone, *n* children gather clockwise in a circle, numbered 0 through *n-1*. Child 0, known as the caller, whispers a message to Child 1 on the right. The message continually passes to Child 2, Child 3... Child *n-1*, and then back to Child 0. As children do not enunciate their words properly (and some may be naughty), the caller's original message is often garbled up into something completely different. 


# The Game of Telephone, with MPI

We want you to implement telephone in MPI. Start with n MPI ranks and a distinguished string. MPI rank 0 sends the string to MPI rank 1, MPI rank 1 then sends string to MPI Rank 2... so on and so forth until the string goes back to MPI rank 0. Before each rank transmits its string, a function *garble()* is applied to the string. We have provided a skeleton code **telephone.c** to be used in the following way:

Let's start with an example:
```
mpirun -n 3 ./telephone "Testing"
```

The first argument (3) denotes the number of MPI ranks and the second argument ("Testing") denotes the string to be passed around. The output should be the following:


```
MPI rank 0 starting message: "Testing"
MPI rank 1 recieved message: "Teshing"
MPI rank 2 recieved message: "Meshing"
MPI rank 0 recieved message: "Mashing"
```

There are three MPI ranks. Rank 0 starts with the original string, and each time the message moves to a different rank it changes due to *garble()*. You should not need any other MPI functions besides *MPI_Send()* and *MPI_Recv()*

More generally, telephone should be used in the following way:

```
mpirun -n NumRanks ./telephone String
```

We have provided a makefile and a submission script. For submission, please take a look at **telephone.c** and finish the telephone implementation.

# Concluding Thoughts
If you are still having trouble, please see the basic Ping Pong program in this [link](http://mpitutorial.com/tutorials/mpi-send-and-receive/). It is quite similar to our game of Telephone, so reading it will provide a few hints and tips!
