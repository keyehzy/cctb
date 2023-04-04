#pragma once

#include <cmath>
#include <cstdlib>
#include <cstring>

template <typename T>
class NBuffer {
 public:
  NBuffer(int size) : size_(size) {
    data_ = static_cast<T*>(std::aligned_alloc(32, size * sizeof(T)));
    std::memset(data_, 0, size * sizeof(T));
  }

  ~NBuffer() { std::free(data_); }

  int size() const { return size_; }

  T* data() const { return data_; }

  T& operator[](int index) { return data_[index]; }

  const T& operator[](int index) const { return data_[index]; }

  NBuffer clone() {
    NBuffer<T> new_buffer(size_);
    std::copy(data_, data_ + size_, new_buffer.data_);
    return new_buffer;
  }

  void swap(NBuffer& other) {
    std::swap(size_, other.size_);
    std::swap(data_, other.data_);
  }

  void fill(T value) { std::fill(data_, data_ + size_, value); }

 private:
  int size_;
  T* data_;
};

template <typename T>
class NumericArray {
 public:
  NumericArray(int size);
  NumericArray(int start, int size, int stride);

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

  NumericArray clone() { return NumericArray(start_, size_, stride_, buffer_.clone()); }

  void swap(NumericArray& other) {
    std::swap(start_, other.start_);
    std::swap(size_, other.size_);
    std::swap(stride_, other.stride_);
    buffer_.swap(other.buffer_);
  }

  NumericArray reshape(int start, int size, int stride) {
    return NumericArray(start, size, stride, buffer_);
  }

  NumericArray flatten() { return NumericArray(0, size_ * stride_, 1, buffer_); }

  T dot(const NumericArray& other) const;
  T norm() const;
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
bool operator==(const NumericArray<T>& lhs, const NumericArray<T>& rhs) {
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
