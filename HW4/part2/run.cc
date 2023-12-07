#include <mpi.h>
#include <cstdio>
void construct_matrices(int *n_ptr, int *m_ptr, int *l_ptr,
                        int **a_mat_ptr, int **b_mat_ptr)
{
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
    {
        scanf("%d %d %d", n_ptr, m_ptr, l_ptr);
        *a_mat_ptr = (int *)calloc((*n_ptr) * (*m_ptr), sizeof(int));
        *b_mat_ptr = (int *)calloc((*m_ptr) * (*l_ptr), sizeof(int));
        for (int i = 0; i < *n_ptr; ++i)
            for (int j = 0; j < *m_ptr; ++j)
                scanf("%d", (*a_mat_ptr) + i * (*m_ptr) + j);
        for (int i = 0; i < *m_ptr; ++i)
            for (int j = 0; j < *l_ptr; ++j)
                scanf("%d", (*b_mat_ptr) + i * (*l_ptr) + j);
    }
}

void matrix_multiply(const int n, const int m, const int l,
                     const int *a_mat, const int *b_mat)
{
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int rows, offset;
    if (world_rank == 0)
    {
        int *c = (int *)calloc(n * l, sizeof(int));
        int averow = n / (world_size - 1);
        int redun = n % (world_size - 1);
        offset = 0;
        for (int rank = 1; rank < world_size; rank++)
        {
            if (rank <= redun)
            {
                rows = averow + 1;
            }
            else
            {
                rows = averow;
            }
            MPI_Send(&n, 1, MPI_INT, rank, 777, MPI_COMM_WORLD);
            MPI_Send(&m, 1, MPI_INT, rank, 777, MPI_COMM_WORLD);
            MPI_Send(&l, 1, MPI_INT, rank, 777, MPI_COMM_WORLD);
            MPI_Send(&offset, 1, MPI_INT, rank, 777, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, rank, 777, MPI_COMM_WORLD);
            MPI_Send(&a_mat[offset * m], rows * m, MPI_INT, rank, 777, MPI_COMM_WORLD);
            MPI_Send(&b_mat[0], m * l, MPI_INT, rank, 777, MPI_COMM_WORLD);
            offset += rows;
        }
        for (int rank = 1; rank < world_size; rank++)
        {
            MPI_Recv(&offset, 1, MPI_INT, rank, 111, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            MPI_Recv(&rows, 1, MPI_INT, rank, 111, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            MPI_Recv(&c[offset * l], rows * l, MPI_INT, rank, 111, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < l; j++)
            {
                printf("%d", c[i * l + j]);
                printf(" ");
            }
            printf("\n");
        }
        free(c);
    }
    else
    {
        int Recive_N, Recive_M, Recive_L;
        MPI_Recv(&Recive_N, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(&Recive_M, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(&Recive_L, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        int *a = (int *)calloc(Recive_N * Recive_M, sizeof(int));
        int *b = (int *)calloc(Recive_M * Recive_L, sizeof(int));
        int *c = (int *)calloc(Recive_N * Recive_L, sizeof(int));
        MPI_Recv(&offset, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(&rows, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(&a[0], rows * Recive_M, MPI_INT, 0, 777, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(&b[0], Recive_M * Recive_L, MPI_INT, 0, 777, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        for (int k = 0; k < Recive_L; k++)
        {
            for (int i = 0; i < rows; i++)
            {
                c[i * Recive_L + k] = 0;
                for (int j = 0; j < Recive_M; j++)
                {
                    c[i * Recive_L + k] += a[i * Recive_M + j] * b[j * Recive_L + k];
                }
            }
        }
        MPI_Send(&offset, 1, MPI_INT, 0, 111, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, 0, 111, MPI_COMM_WORLD);
        MPI_Send(&c[0], rows * Recive_L, MPI_INT, 0, 111, MPI_COMM_WORLD);
        free(a);
        free(b);
        free(c);
    }
}

void destruct_matrices(int *a_mat, int *b_mat)
{
    int word_size, word_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &word_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &word_rank);
    if (word_rank == 0)
    {
        free(a_mat);
        free(b_mat);
    }
}