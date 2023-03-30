#pragma once

#include <cmath>
#include <cstdlib>

#define ALWAYS_INLINE inline __attribute__((always_inline))

template <typename T>
class NBuffer {
 public:
  NBuffer(int size)
      : size_(size), data_(static_cast<T*>(std::aligned_alloc(32, size * sizeof(T)))) {}

  NBuffer(const NBuffer& other) = delete;
  NBuffer& operator=(const NBuffer& other) = delete;

  ~NBuffer() { std::free(data_); }

  int size() const { return size_; }

  T* data() const { return data_; }

  T& operator[](int index) { return data_[index]; }

  const T& operator[](int index) const { return data_[index]; }

 private:
  int size_;
  T* data_;
};

template <typename T>
class NumericArray {
 public:
  NumericArray(int size);
  NumericArray(int start, int size, int stride);
  NumericArray(const NumericArray& other) = delete;
  NumericArray& operator=(const NumericArray& other) = delete;

  ~NumericArray() {}

  T& operator[](int index) { return buffer_[start_ + index * stride_]; }

  const T& operator[](int index) const { return buffer_[start_ + index * stride_]; }

  int start() const { return start_; }

  int size() const { return size_; }

  int stride() const { return stride_; }

  NBuffer<T>& buffer() { return buffer_; }

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
ALWAYS_INLINE bool operator==(const NumericArray<T>& lhs, const NumericArray<T>& rhs) {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (int i = 0; i < lhs.size(); ++i) {
    if (lhs[i] != rhs[i]) {
      return false;
    }
  }
  return true;
}

template <typename T>
ALWAYS_INLINE void axpy(T alpha, const NumericArray<T>& x, NumericArray<T>& y) {
  for (int i = 0; i < x.size(); ++i) {
    y[i] += alpha * x[i];
  }
}
