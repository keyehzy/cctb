#pragma once

class RefCounted {
 public:
  RefCounted() : cnt_(0) {}
  virtual ~RefCounted() {}

  void ref() const { ++cnt_; }

  void unref() const {
    if (--cnt_ == 0) {
      delete this;
    }
  }

 private:
  mutable long cnt_;
};
