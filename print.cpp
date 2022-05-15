#include "stdio.h"
#include "print.hpp"
#include <mpi.h>

void print(double *mat, int n, int m, int l, int k, bool margin) // outputs no more than l x k matrixs
{
    if(margin) printf("\n");
    if (l > n)
        l = n;
    if (k > m)
        k = m;
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < k; j++)
            printf("% 8.2E", mat[i * m + j]);
        printf("\n");
    }
    if(margin) printf("\n\n");
}


void print_mpi(double *mat_part, double *buff,int n,int m,int size,int rank){ // works only with square matrices; outputs no more than m x m matrixs

	MPI_Status status;

	m = (n < m) ? n : m;

	if (rank == 0) printf("\n");
	for (int i = 0; i < m; i++)
	{
		if (rank == 0)
		{
			if (rank == i%size)
			{
				for (int j = 0; j < m; j++){
					printf("% 8.2E ", mat_part[i/size * n + j]);
				}
				printf("\n");
			}
			else
			{
				MPI_Recv(buff, m, MPI_DOUBLE, i%size, 0, MPI_COMM_WORLD, &status);
				for (int j = 0; j < m; j++)
					printf("% 8.2E ", buff[j]);
				printf("\n");
			}
		}
		else if (rank == i%size)
		{
			for (int j = 0; j < m; j++)
				buff[j] = mat_part[i/size * n + j];
			MPI_Send(buff, m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
	if (rank == 0) printf("\n");
}


void save_mpi(double *mat_gather,double *mat_part, double *buff,int n,int size,int rank){
    MPI_Status status;

    for (int i = 0; i < n; i++)
	{
		if (rank == 0)
		{
			if (rank == i%size)
			{
				for (int j = 0; j < n; j++)
					mat_gather[i*n + j] = mat_part[i/size * n + j];
			}
			else
			{
				MPI_Recv(buff, n, MPI_DOUBLE, i%size, 0, MPI_COMM_WORLD, &status);
				for (int j = 0; j < n; j++)
					mat_gather[i*n + j] = buff[j];
			}
		}
		else if (rank == i%size)
		{
			for (int j = 0; j < n; j++)
				buff[j] = mat_part[i/size * n + j];
			MPI_Send(buff, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
}