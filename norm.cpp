#include "norm.hpp"
#include <math.h>
#include <iostream>
using namespace std;

double norm(double* mat, int n, int m) {
    double max = 0, sum = 0;
    for (int j = 0; j < m; ++j) {
        sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += fabs(mat[i * m + j]);
        }
        if (sum > max) max = sum;
    }
    return max;
}

double E_norm(double* mat, int n, int m) { // the Euclidean norm of mat - E 
    double ans = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (i == j)
            {
                ans += (mat[i*m + j]-1)*(mat[i*m + j]-1);
            }
            else
            {
                ans += mat[i*m + j]*mat[i*m + j];
            }
        }
    }
    return sqrt(ans);
}

void substract(double* mat1, double* mat2, double* res, int n, int m) {
    for (int i = 0; i < n * m; ++i) res[i] = mat1[i] - mat2[i];
}

void mult(double* mat1, double* mat2, double* res, int n, int m, int k) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            double c = 0;
            for (int r = 0; r < m; r++) {
                c += mat1[i * m + r] * mat2[r * k + j];
            }
            res[i * k + j] = c;
        }
    }
}