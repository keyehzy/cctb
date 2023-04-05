#pragma once

#include "LinearAlgebra/NumericArray.h"

template <typename T>
class Matrix final : public NumericArray<T> {
 public:
  Matrix(size_t rows, size_t cols) : NumericArray<T>(0, cols, rows) {}

  T& operator()(size_t i, size_t j) { return this->buffer_[i * cols() + j]; }
  T operator()(size_t i, size_t j) const { return this->buffer_[i * cols() + j]; }

  size_t rows() const { return this->stride_; }
  size_t cols() const { return this->size_; }
};
