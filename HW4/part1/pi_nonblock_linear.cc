#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: MPI init
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    long long int per_toss = tosses / world_size;
    unsigned int seed = (unsigned)world_rank * 114;
    srand(seed);
    double distance;
    long long int local_inside = 0;
    long long int total_inside = 0;
    for (long long int i = 0; i < per_toss; i++)
    {
        double x, y;
        x = (double)rand_r(&seed) / RAND_MAX;
        y = (double)rand_r(&seed) / RAND_MAX;

        distance = x * x + y * y;
        if (distance < 1.0)
        {
            local_inside += 1;
        }
    }
    if (world_rank > 0)
    {
        // TODO: MPI workers
        MPI_Request request;
        MPI_Isend( &local_inside ,1 , MPI_LONG_LONG_INT , 0 , 777 , MPI_COMM_WORLD , &request);
    }
    else if (world_rank == 0)
    {
        // TODO: non-blocking MPI communication.
        // Use MPI_Irecv, MPI_Wait or MPI_Waitall.
        MPI_Request requests[world_size-1];
        MPI_Status status[world_size - 1];
        long long int buffer[world_size ];
        total_inside = local_inside;
        for (int i = 1 ; i < world_size;i++){
            MPI_Irecv( &buffer[i] , 1 , MPI_LONG_LONG_INT , i , 777 , MPI_COMM_WORLD , &requests[i-1]);
        
        }

        MPI_Waitall( world_size-1 , requests , MPI_STATUSES_IGNORE);
        for (int i = 1 ; i < world_size;i++){
            total_inside+=buffer[i];
        } 
    }

    if (world_rank == 0)
    {
        // TODO: PI result
        pi_result = 4 * (double)total_inside / (double)tosses;
        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
