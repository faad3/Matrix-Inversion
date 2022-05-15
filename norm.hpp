#ifndef NORM_HPP
#define NORM_HPP

double norm(double* mat, int n, int m);
double E_norm(double* mat, int n, int m);
void substract(double* mat1, double* mat2, double* res, int n, int m);
void mult(double* mat1, double* mat2, double* res, int n, int m, int k);

#endif