#ifndef PRINT_HPP
#define PRINT_HPP

void print(double* mat, int n, int m, int l,int k,bool margin = true);

void print_mpi(double *mat_part, double *buff, int n,int m,int size,int rank);

void save_mpi(double *mat_gather,double *mat_part, double *buff,int n,int size,int rank);

#endif