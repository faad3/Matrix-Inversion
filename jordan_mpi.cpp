#include <iostream>
#include <math.h>
#include "print.hpp"
#include <mpi.h>
#include "jordan_mpi.hpp"

#include <iostream>

using namespace std;

#define EPS 1e-10

void swap_columns_mpi(double *mat_part, int i, int j, int n,int size,int rank)
{
    if(i == j) return;
    int rows;
    if (rank + 1 > n % size)
        rows = n / size;
    else
        rows = n / size + 1;

    double tmp = 0;
    for (int k = 0; k < rows; k++)
    {
        tmp = mat_part[k * n + j];
        mat_part[k * n + j] = mat_part[k * n + i];
        mat_part[k * n + i] = tmp;
    }
}

void swap_rows_mpi(double *mat_part, int i, int j, int n,double *buff,int size,int rank)
{
    MPI_Status status;

    if(i == j) return;
    if(rank == i%size && rank == j%size){
        i = (i-rank)/size;
        j = (j-rank)/size;
        double tmp = 0;
        for (int k = 0; k < n; k++)
        {
            tmp = mat_part[i * n + k];
            mat_part[i * n + k] = mat_part[j * n + k];
            mat_part[j * n + k] = tmp;
        }
    } else if(rank == i%size){
        for(int k = 0; k < n; k++)
            buff[k] = mat_part[n * (i-i%size)/size + k];
        MPI_Recv(&mat_part[n * (i-i%size)/size],n,MPI_DOUBLE,j%size,0,MPI_COMM_WORLD,&status);
        MPI_Send(buff,n,MPI_DOUBLE,j%size,0,MPI_COMM_WORLD);
    } else if(rank == j%size){
        MPI_Send(&mat_part[n * (j-j%size)/size],n,MPI_DOUBLE,i%size,0,MPI_COMM_WORLD);
        MPI_Recv(&mat_part[n * (j-j%size)/size],n,MPI_DOUBLE,i%size,0,MPI_COMM_WORLD,&status);
    }
}   

void find_max_mpi(double &max, int &max_i, int &max_j, double *mat_part, int s, int n, int size, int rank)
{
    int rows;
    if (rank + 1 > n % size)
        rows = n / size;
    else
        rows = n / size + 1;

    int start;
    if (s == 0)
        start = 0;
    else if (rank <= (s - 1) % size)
        start = (s - 1) / size + 1;
    else if (rank > (s - 1) % size)
        start = (s - 1) / size;


    if(rank == (n-1)%size){ 
        max = mat_part[rows*n-1];
        max_i = n-1;
        max_j = n-1;
    }
    MPI_Bcast(&max,1,MPI_DOUBLE,(n-1)%size,MPI_COMM_WORLD);
    MPI_Bcast(&max_i,1,MPI_INT,(n-1)%size,MPI_COMM_WORLD);
    MPI_Bcast(&max_j,1,MPI_INT,(n-1)%size,MPI_COMM_WORLD);

    for (int i = start; i < rows; i++)
    {
        for (int j = s; j < n; j++)
        {
            if (mat_part[i * n + j] > max && fabs(max) > EPS && fabs(mat_part[i * n + j]) > EPS)
            {
                max = mat_part[i * n + j];
                max_i = rank + i * size;
                max_j = j;
            }
            else if (fabs(max) < EPS && fabs(mat_part[i * n + j]) > EPS)
            {
                max = mat_part[i * n + j];
                max_i = rank + i * size;
                max_j = j;
            }
        }
    }
    
    
    dbl_2int max_ij = {max,max_i,max_j};
    MPI_Aint int_ex, double_ex,d_lb,i_lb;
    MPI_Type_get_extent(MPI_INT, &i_lb,&int_ex);
    MPI_Type_get_extent(MPI_DOUBLE,&d_lb,&double_ex);
    int blocklengths[3]={1,1,1};
    MPI_Datatype MPI_DOUBLE_2INT;
    MPI_Datatype types[3]={MPI_DOUBLE, MPI_INT,MPI_INT};
    MPI_Aint displacements[3] = { offsetof(dbl_2int, max),
                                  offsetof(dbl_2int, i),
                                  offsetof(dbl_2int, j)};
    MPI_Type_create_struct(3, blocklengths, displacements, types, &MPI_DOUBLE_2INT);
    MPI_Type_commit(&MPI_DOUBLE_2INT);
    MPI_Op MPI_MAX2LOC;
    MPI_Op_create(max_2loc, 1, &MPI_MAX2LOC);

    MPI_Allreduce(&max_ij,&max_ij,1,MPI_DOUBLE_2INT,MPI_MAX2LOC,MPI_COMM_WORLD);
    max = max_ij.max;
    max_i = max_ij.i;
    max_j = max_ij.j;

}



void max_2loc(void *in, void *inout, int *len, MPI_Datatype *type){
    /* ignore type, just trust that it's our dbl_twoindex type */
    dbl_2int *invals    = (dbl_2int*) in;
    dbl_2int *inoutvals = (dbl_2int*) inout;

    for (int i=0; i < *len; i++) {
        if (invals[i].max > inoutvals[i].max) {
            inoutvals[i].max  = invals[i].max;
            inoutvals[i].i = invals[i].i;
            inoutvals[i].j = invals[i].j;
        }
    }

    return;
}

void mult_row_mpi(double *mat_part, int i, double mult, int n,int size,int rank)
{
    if(rank == i%size){
        i = (i-rank)/size;
        for (int j = 0; j < n; ++j)
            mat_part[i * n + j] *= mult;
    }
}

void substr_row(double *mat, int i, int j, double mult, int m)
{
    for (int k = 0; k < m; ++k)
        mat[j * m + k] -= mult * mat[i * m + k];
}

void substr_mpi(double *mat_part,double *inv_mat_part,int n,int s, double *buff,int size,int rank){
    double mult;

    int rows;
    if (rank + 1 > n % size)
        rows = n / size;
    else
        rows = n / size + 1;


    double *mat_buff = buff;
    double *inv_buff = buff+n;
    if(rank == s%size){
        for(int k =0;k<n;k++){
            mat_buff[k] = mat_part[n * (s-rank)/size+k];
            inv_buff[k] = inv_mat_part[n * (s-rank)/size+k];
        }
    } 
    MPI_Bcast(mat_buff,n,MPI_DOUBLE,s%size,MPI_COMM_WORLD);
    MPI_Bcast(inv_buff,n,MPI_DOUBLE,s%size,MPI_COMM_WORLD);
    for(int i = 0; i < rows;i++){
        if(rank + i*size != s){
            mult = mat_part[i*n+s];
            for (int k = 0; k < n; k++){
                mat_part[i * n + k] -= mult * mat_buff[k];
                inv_mat_part[i * n + k] -= mult * inv_buff[k];
            }
        }
    }
}

bool inverse_mpi(double *mat_part, double *inv_mat_part, int *col_swaps, int n,double *buff,int size,int rank)
{
    int max_i = 0;
    int max_j = 0;
    double max;
    double mult;
    for (int s = 0; s < n; s++)
    { 
        find_max_mpi(max, max_i, max_j, mat_part, s,  n,size,rank);
        if (fabs(max) < EPS)
        {
            if(rank == 0) cout <<  "Degenerate matrix!" << endl;
            return 0;
        }

        swap_rows_mpi(mat_part, s, max_i, n,buff,size,rank);
        swap_rows_mpi(inv_mat_part, s, max_i, n,buff,size,rank);

        swap_columns_mpi(mat_part, s, max_j, n,size,rank);
        swap_columns_mpi(inv_mat_part, s, max_j, n,size,rank);
        col_swaps[s] = max_j;

        mult_row_mpi(mat_part, s, 1.0 / max, n,size,rank);
        mult_row_mpi(inv_mat_part, s, 1.0 / max, n,size,rank);

        substr_mpi(mat_part,inv_mat_part,n,s,buff,size,rank);
    }



    for (int s = n - 1; s >= 0; s--)
    {
        if (s != col_swaps[s])
        {
            swap_columns_mpi(inv_mat_part, s, col_swaps[s], n,size,rank);
            swap_columns_mpi(mat_part, s, col_swaps[s], n,size,rank);
        }
    }

     
    

    int rows;
    if (rank + 1 > n % size)
        rows = n / size;
    else
        rows = n / size + 1;

    int t = 0;
    int t_rank[2];
    for(int j = 0; j < n-1; j++){
        t = 0;
        t_rank[1] = rank;
        while ((t < rows-1 ) && (fabs(mat_part[t * n + j] - 1) > EPS)){
            t++;
        }
        if((fabs(mat_part[t * n + j] - 1) > EPS)) t = -1;
        t_rank[0] = t;
        MPI_Allreduce(&t_rank,&t_rank,1,MPI_2INT,MPI_MAXLOC,MPI_COMM_WORLD);

        swap_rows_mpi(mat_part, j, t_rank[1]+t_rank[0]*size, n,buff,size,rank);
        swap_rows_mpi(inv_mat_part, j, t_rank[1]+t_rank[0]*size, n,buff,size,rank);
    }


    return 1;
}