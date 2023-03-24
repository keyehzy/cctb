#pragma once

template <typename T>
struct Matrix {
  int rows;
  int cols;
  T* data;

  Matrix(int rows, int cols) : rows(rows), cols(cols) {
    data = new T[rows * cols];
  }

  ~Matrix() { delete[] data; }

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
};