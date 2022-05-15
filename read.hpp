#include <string>

using namespace std;

#ifndef READ_HPP
#define READ_HPP

double func(int k, int n, int i, int j);

void read_func(double* mat, int n, int k);
void read_file(double* mat, int n, string name);

void read_func_mpi(double* mat_part, int n, int k,int size,int rank);
void read_file_mpi(double* mat_part, int n, string name,int size,int rank);

void E_mat_mpi(double *mat_part,int n,int size,int rank);

#endif