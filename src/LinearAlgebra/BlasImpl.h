#pragma once

#include <complex>

#include "LinearAlgebra/Matrix.h"
#include "LinearAlgebra/NumericArray.h"

void scal(double alpha, NumericArray<double>& x);
void copy(const NumericArray<double>& x, NumericArray<double>& y);
void axpy(double alpha, const NumericArray<double>& x, NumericArray<double>& y);
void dot(const NumericArray<double>& x, const NumericArray<double>& y, double* result);
void nrm2(const NumericArray<double>& x, double* result);
void asum(const NumericArray<double>& x, double* result);
void iamax(const NumericArray<double>& x, int* result);

void scal(std::complex<double> alpha, NumericArray<std::complex<double>>& x);
void copy(const NumericArray<std::complex<double>>& x, NumericArray<std::complex<double>>& y);
void axpy(std::complex<double> alpha, const NumericArray<std::complex<double>>& x,
          NumericArray<std::complex<double>>& y);
void dot(const NumericArray<std::complex<double>>& x, const NumericArray<std::complex<double>>& y,
         std::complex<double>* result);
void asum(const NumericArray<std::complex<double>>& x, std::complex<double>* result);
void iamax(const NumericArray<std::complex<double>>& x, int* result);

void gemv(double alpha, const Matrix<double>& A, const NumericArray<double>& x, double beta,
          NumericArray<double>& y);
void ger(double alpha, const NumericArray<double>& x, const NumericArray<double>& y,
         Matrix<double>& A);
void gemm(double alpha, const Matrix<double>& a, const Matrix<double>& b, double beta,
          Matrix<double>& c);

void gemv(std::complex<double> alpha, const Matrix<std::complex<double>>& A,
          const NumericArray<std::complex<double>>& x, std::complex<double> beta,
          NumericArray<std::complex<double>>& y);
void ger(std::complex<double> alpha, const NumericArray<std::complex<double>>& x,
         const NumericArray<std::complex<double>>& y, Matrix<std::complex<double>>& A);
void gemm(std::complex<double> alpha, const Matrix<std::complex<double>>& a,
          const Matrix<std::complex<double>>& b, std::complex<double> beta,
          Matrix<std::complex<double>>& c);
