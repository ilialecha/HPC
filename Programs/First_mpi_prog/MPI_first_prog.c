#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int message = world_rank + 1;
    int received_message;

    // Point-to-point communication
    if (world_rank == 0) {
        MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&received_message, 1, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&received_message, 1, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&message, 1, MPI_INT, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD);
    }

    printf("Process %d received message %d from process %d\n", world_rank, received_message, (world_rank - 1 + world_size) % world_size);

    // Collective communication
    int sum;
    MPI_Allreduce(&message, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    printf("Process %d calculated the sum of all messages: %d\n", world_rank, sum);

    MPI_Finalize();
    return 0;
}