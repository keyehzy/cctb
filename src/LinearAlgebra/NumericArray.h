#pragma once

#include <cmath>
#include <cstdint>

#include "LinearAlgebra/ContiguousArray.h"

#define ALWAYS_INLINE inline __attribute__((always_inline))

template <typename T>
class NumericArray {
 public:
  NumericArray(int size);
  NumericArray(int start, int size, int stride);
  NumericArray(ContiguousArray<T> data, int start, int size, int stride);
  NumericArray(const NumericArray& other);

  virtual ~NumericArray() {}

  T& operator[](int index) { return contiguous_->at(start_ + index * stride_); }

  const T& operator[](int index) const { return contiguous_->at(start_ + index * stride_); }

  NumericArray& operator=(const NumericArray& other);

  NumericArray into(int start, int size, int stride) const;

  int size() const { return size_; }

  int stride() const { return stride_; }

  T* data() const { return contiguous_->data() + start_; }

  template <typename U>
  friend ALWAYS_INLINE bool operator==(const NumericArray<U>& lhs, const NumericArray<U>& rhs);

 protected:
  int start_;
  int size_;
  int stride_;
  RefPtr<ContiguousArray<T>> contiguous_;
};

template <typename T>
NumericArray<T>::NumericArray(int size)
    : start_(0), size_(size), stride_(1), contiguous_(new ContiguousArray<T>(size, 32)) {}

template <typename T>
NumericArray<T>::NumericArray(int start, int size, int stride)
    : start_(start),
      size_(size),
      stride_(stride),
      contiguous_(new ContiguousArray<T>(size * stride * sizeof(T), 32)) {}

template <typename T>
NumericArray<T>::NumericArray(ContiguousArray<T> data, int start, int size, int stride)
    : start_(start), size_(size), stride_(stride), contiguous_(data) {}

template <typename T>
NumericArray<T>::NumericArray(const NumericArray& other)
    : start_(other.start_),
      size_(other.size_),
      stride_(other.stride_),
      contiguous_(other.contiguous_) {}

template <typename T>
NumericArray<T> NumericArray<T>::into(int start, int size, int stride) const {
  return NumericArray<T>(contiguous_, start, size, stride);
}

template <typename T>
NumericArray<T>& NumericArray<T>::operator=(const NumericArray& other) {
  if (this != &other) {
    contiguous_ = other.contiguous_;
    start_ = other.start_;
    size_ = other.size_;
    stride_ = other.stride_;
  }
  return *this;
}
