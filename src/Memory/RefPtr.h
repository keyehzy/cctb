#pragma once

template <class RefCountedObject>
class RefPtr {
 public:
  RefPtr() : obj_(nullptr) {}

  RefPtr(const RefPtr& other) { add_ref(other.obj_); }

  RefPtr(RefCountedObject* obj) { add_ref(obj); }

  ~RefPtr() { remove_ref(); }

  void operator=(const RefPtr& other) {
    remove_ref();
    add_ref(other.obj_);
  }

  void operator=(RefCountedObject* obj) {
    remove_ref();
    add_ref(obj);
  }

  RefCountedObject* operator->() const { return obj_; }

 private:
  void add_ref(RefCountedObject* new_obj) {
    obj_ = new_obj;
    if (obj_ != nullptr) {
      obj_->ref();
    }
  }

  void remove_ref() {
    if (obj_ != nullptr) {
      obj_->unref();
    }
  }

  RefCountedObject* obj_;
};