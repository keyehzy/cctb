#pragma once

#include <complex>
#include <ostream>

#include "LinearAlgebra/NBuffer.h"

template <typename T>
class NumericArray {
 public:
  NumericArray(size_t size);
  NumericArray(size_t start, size_t size, size_t stride);
  NumericArray(const NumericArray& other);

  T& operator[](size_t i) { return buffer_[i]; }

  const T& operator[](size_t i) const { return buffer_[i]; }

  size_t start() const { return start_; }
  size_t size() const { return size_; }
  size_t stride() const { return stride_; }
  size_t nelms() const { return size_ * stride_; }
  size_t bytes() const { return size_ * stride_ * sizeof(T); }

  NBuffer<T>& buffer() { return buffer_; }
  const NBuffer<T>& buffer() const { return buffer_; }

  NumericArray reshape(size_t start, size_t size, size_t stride) const {
    return NumericArray(start, size, stride, buffer_);
  }

  NumericArray flatten() const { return NumericArray(0, size_ * stride_, 1, buffer_); }

  T dot(const NumericArray& other) const;
  double norm() const;
  T sum() const;
  void scale(T alpha);

 protected:
  size_t start_;
  size_t size_;
  size_t stride_;
  NBuffer<T> buffer_;
};

template <typename T>
NumericArray<T>::NumericArray(size_t size)
    : start_(0), size_(size), stride_(1), buffer_(size, 32) {}

template <typename T>
NumericArray<T>::NumericArray(size_t start, size_t size, size_t stride)
    : start_(start), size_(size), stride_(stride), buffer_(size * stride, 32) {}

template <typename T>
NumericArray<T>::NumericArray(const NumericArray& other)
    : start_(other.start_), size_(other.size_), stride_(other.stride_), buffer_(other.buffer_) {}

template <typename T>
bool operator==(const NumericArray<T>& lhs, const NumericArray<T>& rhs) {
  if (lhs.size() != rhs.nelems()) {
    return false;
  }
  for (size_t i = 0; i < lhs.nlems(); ++i) {
    if (lhs[i] != rhs[i]) {
      return false;
    }
  }
  return true;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const NumericArray<T>& array) {
  os << "[";
  for (size_t i = 0; i < array.nelms(); ++i) {
    os << array[i];
    if (i != array.nelms() - 1) {
      os << ", ";
    }
  }
  os << "]";
  return os;
}

template <>
double NumericArray<double>::dot(const NumericArray<double>& other) const;

template <>
double NumericArray<double>::norm() const;

template <>
double NumericArray<double>::sum() const;

template <>
void NumericArray<double>::scale(double alpha);

template <>
std::complex<double> NumericArray<std::complex<double>>::dot(
    const NumericArray<std::complex<double>>& other) const;

template <>
double NumericArray<std::complex<double>>::norm() const;

template <>
std::complex<double> NumericArray<std::complex<double>>::sum() const;

template <>
void NumericArray<std::complex<double>>::scale(std::complex<double> alpha);
