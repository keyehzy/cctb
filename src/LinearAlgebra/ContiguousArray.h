#include <cstdlib>

#include "Memory/RefCounted.h"
#include "Memory/RefPtr.h"

template <typename T>
class ContiguousArray : public RefCounted {
 public:
  ContiguousArray(int size, int align) : size_(size), align_(align) {
    data_ = static_cast<T*>(std::aligned_alloc(align, size * sizeof(T)));
  }

  ~ContiguousArray() { std::free(data_); }

  T at(int index) const { return data_[index]; }

  T& at(int index) { return data_[index]; }

  T* data() const { return data_; }

  int size() const { return size_; }

  int align() const { return align_; }

 private:
  T* data_;
  int size_;
  int align_;
};