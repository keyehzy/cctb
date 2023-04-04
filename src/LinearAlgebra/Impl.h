#pragma once

#include <complex>

#include "LinearAlgebra/NumericArray.h"

void scal(double alpha, NumericArray<double>& x);
void scal(std::complex<double> alpha, NumericArray<std::complex<double>>& x);
void copy(const NumericArray<double>& x, NumericArray<double>& y);
void copy(const NumericArray<std::complex<double>>& x, NumericArray<std::complex<double>>& y);
void axpy(double alpha, const NumericArray<double>& x, NumericArray<double>& y);
void axpy(std::complex<double> alpha, const NumericArray<std::complex<double>>& x,
          NumericArray<std::complex<double>>& y);
void dot(const NumericArray<double>& x, const NumericArray<double>& y, double* result);
void dot(const NumericArray<std::complex<double>>& x, const NumericArray<std::complex<double>>& y,
         std::complex<double>* result);
void nrm2(const NumericArray<double>& x, double* result);
void asum(const NumericArray<double>& x, double* result);
void asum(const NumericArray<std::complex<double>>& x, std::complex<double>* result);
void iamax(const NumericArray<double>& x, int* result);
void iamax(const NumericArray<std::complex<double>>& x, int* result);
