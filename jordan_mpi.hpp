#include <mpi.h>

#ifndef JORDAN_MPI_HPP
#define JORDAN_MPI_HPP

typedef struct {
    double max;
    int i;
    int j;
} dbl_2int;
void max_2loc(void *in, void *inout, int *len, MPI_Datatype *type);


void swap_columns_mpi(double *mat, int i, int j, int n,int size,int rank);
void swap_rows_mpi(double *mat_part, int i, int j, int n,double *buff,int size,int rank);
void find_max_mpi(double &max, int &max_i, int &max_j, double *mat_part, int s, int n, int size, int rank);
void mult_row_mpi(double *mat_part, int i, double mult, int n,int size,int rank);
void substr_row(double* mat, int i, int j, double mult, int m);
void substr_mpi(double *mat_part,double *inv_mat_part,int n,int s, double *buff,int size,int rank);
bool inverse_mpi(double *mat_part, double *inv_mat_part, int *col_swaps, int n,double *buff,int size,int rank);


#endif 