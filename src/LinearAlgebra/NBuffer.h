#pragma once

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>

template <typename T>
class NBuffer {
 public:
  NBuffer(size_t size, size_t align) : size_(size), align_(align) {
    size_t bytes = size * sizeof(T);
    size_t next_multiple_of_align = ((bytes + align - 1) / align) * align;
    data_ = static_cast<T*>(std::aligned_alloc(align, next_multiple_of_align));
    std::memset(data_, 0, bytes);
  }

  NBuffer(const NBuffer& other) : size_(0), align_(0), data_(nullptr) {
    NBuffer<T> new_buffer(other.size_, other.align_);
    std::swap(size_, new_buffer.size_);
    std::swap(align_, new_buffer.align_);
    std::swap(data_, new_buffer.data_);
  }

  ~NBuffer() { std::free(data_); }

  size_t size() const { return size_; }

  T* data() const { return data_; }

  T& operator[](size_t index) { return data_[index]; }

  const T& operator[](size_t index) const { return data_[index]; }

  NBuffer clone() const {
    NBuffer<T> new_buffer(size_, align_);
    std::copy(data_, data_ + size_, new_buffer.data_);
    return new_buffer;
  }

  void swap(NBuffer& other) {
    std::swap(size_, other.size_);
    std::swap(align_, other.align_);
    std::swap(data_, other.data_);
  }

  void fill(T value) { std::fill(data_, data_ + size_, value); }

 private:
  size_t size_;
  size_t align_;
  T* data_;
};
