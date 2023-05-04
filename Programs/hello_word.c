/* MPI Lab 1, Example Program   */

#include "mpi.h"
#include <stdio.h>

int main(argc, argv)
int argc;
char **argv;
{
    int rank, size;
    // Initialize MPI 
    MPI_Init(&argc,&argv);
    // Get the rank of this process (process ID)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /*
    * Print a hello world message from each process
    * The message should include the rank of the process
    * and the total number of processes.
    */
    printf("Hello world! I am %d of %d\n",rank,size);
    // Terminate MPI
    MPI_Finalize();

    return 0;
}