#pragma once

#include "LinearAlgebra/NumericArray.h"

template <typename T>
class Matrix final : public NumericArray<T> {
 public:
  Matrix(int rows, int cols) : NumericArray<T>(0, cols, rows) {}

  T& operator()(int i, int j) { return this->buffer_[i * cols() + j]; }
  T operator()(int i, int j) const { return this->buffer_[i * cols() + j]; }

  int rows() const { return this->stride_; }
  int cols() const { return this->size_; }
};

void sgemv(float alpha, const Matrix<float>& A, const NumericArray<float>& x, float beta,
           NumericArray<float>& y);

void sger(float alpha, const NumericArray<float>& x, const NumericArray<float>& y,
          Matrix<float>& A);

void sgemm(float alpha, const Matrix<float>& a, const Matrix<float>& b, float beta,
           Matrix<float>& c);

void dgemv(double alpha, const Matrix<double>& A, const NumericArray<double>& x, double beta,
           NumericArray<double>& y);

void dger(double alpha, const NumericArray<double>& x, const NumericArray<double>& y,
          Matrix<double>& A);

void dgemm(double alpha, const Matrix<double>& a, const Matrix<double>& b, double beta,
           Matrix<double>& c);
