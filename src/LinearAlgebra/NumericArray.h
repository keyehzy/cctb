#pragma once

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>

#define ALWAYS_INLINE inline __attribute__((always_inline))

template <typename T>
class NumericArray {
 public:
  NumericArray(int size);
  NumericArray(int start, int size, int stride);
  NumericArray(T* data, int start, int size, int stride);
  NumericArray(const NumericArray& other);

  ~NumericArray() { std::free(data_); }

  T& operator[](int index) { return data_[start_ + index * stride_]; }

  const T& operator[](int index) const { return data_[start_ + index * stride_]; }

  NumericArray& operator=(const NumericArray& other);

  NumericArray into(int start, int size, int stride) const;

  friend ALWAYS_INLINE bool operator==(const NumericArray<T>& lhs, const NumericArray<T>& rhs);

 protected:
  T* data_;
  int start_;
  int size_;
  int stride_;
};

template <typename T>
NumericArray<T>::NumericArray(int size) : start_(0), size_(size), stride_(1) {
  data_ = (T*)std::aligned_alloc(32, size_ * stride_ * sizeof(T));
}

template <typename T>
NumericArray<T>::NumericArray(int start, int size, int stride)
    : start_(start), size_(size), stride_(stride) {
  data_ = (T*)std::aligned_alloc(32, size_ * stride_ * sizeof(T));
}

template <typename T>
NumericArray<T>::NumericArray(T* data, int start, int size, int stride)
    : data_(data), start_(start), size_(size), stride_(stride) {}

template <typename T>
NumericArray<T>::NumericArray(const NumericArray& other)
    : data_(other.data_), start_(other.start_), size_(other.size_), stride_(other.stride_) {}

template <typename T>
NumericArray<T> NumericArray<T>::into(int start, int size, int stride) const {
  return NumericArray<T>(data_, start, size, stride);
}

template <typename T>
NumericArray<T>& NumericArray<T>::operator=(const NumericArray& other) {
  if (this != &other) {
    data_ = other.data_;
    start_ = other.start_;
    size_ = other.size_;
    stride_ = other.stride_;
  }
  return *this;
}
