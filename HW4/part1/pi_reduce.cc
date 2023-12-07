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
    // TODO: use MPI_Reduce
    MPI_Reduce( &local_inside , &total_inside , 1 , MPI_LONG_LONG_INT , MPI_SUM , 0 , MPI_COMM_WORLD);
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
