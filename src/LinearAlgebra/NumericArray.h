#pragma once

#include <ostream>

#include "LinearAlgebra/NBuffer.h"

template <typename T>
class NumericArray {
 public:
  NumericArray(int size);
  NumericArray(int start, int size, int stride);
  NumericArray(int start, int size, int stride, NBuffer<T> const& buffer);
  NumericArray(NumericArray const& other);
  NumericArray(NumericArray&& other);

  ~NumericArray() {}

  T& operator[](int i) { return buffer_[i]; }

  const T& operator[](int i) const { return buffer_[i]; }

  int start() const { return start_; }
  int size() const { return size_; }
  int stride() const { return stride_; }
  int nelms() const { return size_ * stride_; }
  int bytes() const { return size_ * stride_ * sizeof(T); }

  NBuffer<T>& buffer() { return buffer_; }
  const NBuffer<T>& buffer() const { return buffer_; }

  NumericArray reshape(int start, int size, int stride) const {
    return NumericArray(start, size, stride, buffer_);
  }

  NumericArray flatten() const { return NumericArray(0, size_ * stride_, 1, buffer_); }

  T dot(const NumericArray& other) const;
  double norm() const;
  T sum() const;
  void scale(T alpha);

 protected:
  int start_;
  int size_;
  int stride_;
  NBuffer<T> buffer_;
};

template <typename T>
NumericArray<T>::NumericArray(int size) : start_(0), size_(size), stride_(1), buffer_(size) {}

template <typename T>
NumericArray<T>::NumericArray(int start, int size, int stride)
    : start_(start), size_(size), stride_(stride), buffer_(size * stride) {}

template <typename T>
NumericArray<T>::NumericArray(int start, int size, int stride, NBuffer<T> const& buffer)
    : start_(start), size_(size), stride_(stride), buffer_(buffer) {}

template <typename T>
NumericArray<T>::NumericArray(NumericArray const& other)
    : start_(other.start_), size_(other.size_), stride_(other.stride_), buffer_(other.buffer_) {}

template <typename T>
NumericArray<T>::NumericArray(NumericArray&& other)
    : start_(other.start_),
      size_(other.size_),
      stride_(other.stride_),
      buffer_(std::move(other.buffer_)) {}

template <typename T>
bool operator==(const NumericArray<T>& lhs, const NumericArray<T>& rhs) {
  if (lhs.size() != rhs.nelems()) {
    return false;
  }
  for (int i = 0; i < lhs.nlems(); ++i) {
    if (lhs[i] != rhs[i]) {
      return false;
    }
  }
  return true;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const NumericArray<T>& array) {
  os << "[";
  for (int i = 0; i < array.nelms(); ++i) {
    os << array[i];
    if (i != array.nelms() - 1) {
      os << ", ";
    }
  }
  os << "]";
  return os;
}
