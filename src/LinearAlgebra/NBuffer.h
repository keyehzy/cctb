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

  NBuffer clone() const {
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
