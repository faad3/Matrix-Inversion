#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include "read.hpp"
using namespace std;

double func(int k, int n, int i, int j) 
{
    double ans = 0;
    switch(k){
        case 1:
            ans = n-max(i,j) + 1;
            break;
        case 2:
            ans = max(i,j);
            break;
        case 3:
            ans = fabs(i-j);
            break;
        case 4:
            ans = (double)1/(i+j-1);
            break;
    }
    return ans;
}


void read_func(double* mat, int n, int k) 
{
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            mat[n * i + j] = func(k, n, i+1, j+1);
        }
    }
}

void read_file(double* mat, int n, string name) 
{
    ifstream in;
    in.open(name);
    if (!in.is_open()) 
    {
        cout << "error opening the file" << endl;
        return;
    }
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            in >> mat[i * n + j];
        }
    }
    in.close();
}


void read_func_mpi(double* mat_part, int n, int k,int size,int rank){

    int rows;
    if(rank + 1 > n%size) rows = n/size; 
    else rows = n/size + 1;

    for (int i = 0; i < rows; i++)
		for (int j = 0; j < n; j++)
            mat_part[i*n + j] = func(k, n, rank+size*i + 1, j + 1);
}

void read_file_mpi(double* mat_part, int n, string name,int size,int rank){

    int rows; // number of rows of matrix_part
    if(rank + 1 > n%size) rows = n/size; 
    else rows = n/size + 1;

    ifstream in;
    in.open(name);
    if (!in.is_open()) 
    {
        cout << "error opening the file" << endl;
        return;
    }
    for (int i = 0; i < rows; i++) 
    {
        for (int j = 0; j < n; j ++) 
        {
            in.ignore(rank + i*size);
            in >> mat_part[i*n + j];
        }
    }
    in.close();
}


void E_mat_mpi(double *mat_part,int n,int size,int rank){

    int rows;
    if(rank + 1 > n%size) rows = n/size; 
    else rows = n/size + 1;
    for (int i = 0; i < rows; i++)
		for (int j = 0; j < n; j++)
            mat_part[i*n+j] = ( j == rank+size*i) ? 1 : 0;
}
