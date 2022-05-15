#include "print.hpp"
#include "read.hpp"
#include "jordan_mpi.hpp"
#include "norm.hpp"
#include <time.h>
#include <cmath>
#include <iostream>
#include <string>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(&argc,&argv);

    MPI_Status status;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n, m, k;
    if (argc < 4){
        if (rank == 0) cout << "incorrect number of arguments" << endl;
        MPI_Finalize();
        return -1;
    }

    n = stoi(argv[1]);
    m = stoi(argv[2]);
    k = stoi(argv[3]);
    if (k < 0 || k > 4)
    {
        if (rank == 0) cout << "incorrect argument" << endl;
        MPI_Finalize();
        return -1;
    }
        
    // INPUT (the matrix is divided between processes by rows)

    double *mat_check,*inv_mat_gather; // for validation
    if (rank == 0)  
    {
        mat_check = new double[n * n];
        inv_mat_gather = new double[n * n];
        if (argc == 5 && k == 0)
        {
            string name = argv[4];
            read_file(mat_check, n, name);
        }
        else if (argc == 4 || k == 0)
        {
            read_func(mat_check, n, k);
        }
    }


    int rows; // number of rows matrix_part
    if(rank + 1 > n%size) rows = n/size; 
    else rows = n/size + 1;
    
    double *buff = new double[2*n];
    int* col_swaps = new int[n];
    for(int i =0;i<n;i++)
        col_swaps[i] = i;

    double *mat_part = new double[rows * n];
    double *inv_mat_part = new double[rows * n];
    E_mat_mpi(inv_mat_part,n,size,rank);
    if (argc == 5 && k == 0)
    {
        string name = argv[4];
        read_file_mpi(mat_part, n, name, size, rank);
    }
    else if (argc == 4)
    {
        read_func_mpi(mat_part, n, k, size, rank);
    }
    
    if(rank == 0) cout << "Input matrix: ";
    print_mpi(mat_part,buff,n,m,size,rank);


    // MATRIX INVERSION
    

    MPI_Barrier(MPI_COMM_WORLD);  
    double time = MPI_Wtime();

    inverse_mpi(mat_part,inv_mat_part, col_swaps, n, buff, size,rank);

    MPI_Barrier(MPI_COMM_WORLD);
	time = MPI_Wtime() - time;


    // VALIDATION
    if(rank == 0) cout << "Inverse Matrix:" << endl; 
    print_mpi(inv_mat_part,buff,n,m,size,rank);
    save_mpi(inv_mat_gather,inv_mat_part,buff,n,size,rank);
    double* mat_mult;
    if(rank == 0){
        mat_mult = new double[n*n];
        mult(mat_check,inv_mat_gather,mat_mult,n,n,n);
        cout << "Multipy result:" << endl;
        print(mat_mult, n, n, m,m);   
        cout << "Residual: " << E_norm(mat_mult,n,n) << endl;  
        cout << "Time :" << time << " sec" << endl << endl;
    }


    // CLEAR MEMORY

    if(rank == 0) {delete[] mat_check; delete[] inv_mat_gather; delete[] mat_mult;}
    delete[] mat_part;delete[] inv_mat_part; delete[] buff;

    MPI_Finalize();
}