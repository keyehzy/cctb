#include "LinearAlgebra/Matrix.h"

void gemv(double alpha, const Matrix<double>& A, const NumericArray<double>& x, double beta,
          NumericArray<double>& y) {
  for (int i = 0; i < A.rows(); ++i) {
    double result = 0;
    for (int j = 0; j < A.cols(); ++j) {
      result += A(i, j) * x[j];
    }
    y[i] = alpha * result + beta * y[i];
  }
}

void ger(double alpha, const NumericArray<double>& x, const NumericArray<double>& y,
         Matrix<double>& A) {
  for (int i = 0; i < A.rows(); ++i) {
    for (int j = 0; j < A.cols(); ++j) {
      A(i, j) += alpha * x[i] * y[j];
    }
  }
}

void gemm(double alpha, const Matrix<double>& a, const Matrix<double>& b, double beta,
          Matrix<double>& c) {
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < b.cols(); ++j) {
      double result = 0;
      for (int k = 0; k < a.cols(); ++k) {
        result += a(i, k) * b(k, j);
      }
      c(i, j) = alpha * result + beta * c(i, j);
    }
  }
}
