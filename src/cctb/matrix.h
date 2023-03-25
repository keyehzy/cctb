#pragma once

#include <iostream>
#include <vector>

template <typename T>
struct Matrix {
  int rows;
  int cols;
  std::vector<T> data;

  Matrix(int rows, int cols) : rows(rows), cols(cols), data(rows * cols) {}

  ~Matrix() {}

  T& operator()(int row, int col) { return data[row * cols + col]; }

  T operator()(int row, int col) const { return data[row * cols + col]; }

  Matrix operator*(T scalar) {
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        result(i, j) = (*this)(i, j) * scalar;
      }
    }
    return result;
  }

  Matrix Print() {
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        std::cout << (*this)(i, j) << " ";
      }
      std::cout << std::endl;
    }
    return *this;
  }
};