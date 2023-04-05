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
