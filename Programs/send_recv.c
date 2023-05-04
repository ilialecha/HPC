/* MPI Lab 2a, Example Program   */
// This program is a simple example of sending and receiving messages
// between processes.  The program is designed to be run with 2 processes.
// The program will send a message from process 0 to process 1.  Process 1
// will then print the message it received.

#include "mpi.h"
#include <stdio.h>

int main(argc, argv)
int argc;
char **argv;
{
    /*
    * Declare variables
    * n is the number of messages to send
    * rank is the process ID
    * size is the number of processes
    * to is the process to send to
    * from is the process to receive from
    * tagno is the message tag number
    * status is the status of the receive
    * message is the message to send
    */
    int rank, size, n, to, from, tagno;
    // Declare status variable for receive operation.
    MPI_Status status;
    // n is the number of messages to send (from process 0)
    n = -1;
    // Initialize MPI
    MPI_Init(&argc,&argv);
    // Get the rank of this process (process ID)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Set the process to send to
    // Process 0 sends to process 1
    to = rank + 1;
    // Check if process is the last process
    // If so, send to process 0 instead
    if (rank == size -1) to = 0;
    // Set the process to receive from
    // Process 0 receives from process 1
    from = rank - 1;
    // tagno is the message tag number (unique identifier)
    tagno = 201;
    // print a message from each process.
    printf("Process  %d of %d is alive\n",rank,size);
    // Check if process is process 0 (sender) or process 1 (receiver)
    if (rank == 0){
        // Process 0 sends a message to process 1
        from = size - 1;

        printf("Please enter a positive integer\n");
        // Get the number of messages to send
        scanf("%d",&n);
        // Print the number of messages to send
        printf("n = %d\n",n);
        // Send the number of messages to process 1
        MPI_Send(&n,1,MPI_INT,to,tagno,MPI_COMM_WORLD);
    }
    // While the number of messages is greater than 0
    while (1){
        // MPI_ANY_SOURCE is a wildcard for the source process
        from = MPI_ANY_SOURCE;
        // Receive a message from any process
        MPI_Recv(&n,1,MPI_INT,from,tagno,MPI_COMM_WORLD, &status); 
        // Print the message received
        printf ("Rank %d received %d\n",rank, n);
        // Check if process is process 0 (sender) or process 1 (receiver)
        // If process is process 0, decrement the number of messages
        // to send and increment the message tag number
        if (rank == 0) {n--;tagno++;}
        // Send the number of messages to process 1
        MPI_Send(&n,1,MPI_INT,to,tagno,MPI_COMM_WORLD);
        // Check if process is process 1 (receiver)
        // If process is process 1, decrement the number of messages
        // to send and increment the message tag number
        if (rank != 0) {n--;tagno++;}
        // Check if the number of messages is less than 0
        if (n<0){
            // Terminate MPI
            MPI_Finalize();
            return 0;
        }
    }
}