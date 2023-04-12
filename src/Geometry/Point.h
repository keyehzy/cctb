#pragma once

#include <array>
#include <iostream>

#include "Geometry/Vector.h"

template <size_t D>
class Vector;
template <size_t D>
class Point;

template <size_t D>
bool operator==(const Point<D>& lhs, const Point<D>& rhs);

template <size_t D>
class Point {
 public:
  Point() : data_({}) {}
  Point(double x) : data_({x}) {}
  Point(double x, double y) : data_({x, y}) {}
  Point(double x, double y, double z) : data_({x, y, z}) {}
  Point(const std::array<double, D>& a) : data_(a) {}
  Point(const Point& other) : data_(other.data_) {}

  int size() const { return D; }

  Point translated(const Vector<D>& v) const {
    Point result;
    for (size_t i = 0; i < D; i++) {
      result.data_[i] = data_[i] + v[i];
    }
    return result;
  }

  Vector<D> as_vector_from_origin() const { return Vector<D>(data_); }

  double distance_to(const Point& other) const { return Vector(*this, other).norm(); }

  double& operator[](size_t i) { return data_[i]; }
  double operator[](size_t i) const { return data_[i]; }

  Point& operator=(const Point& other) {
    data_ = other.data_;
    return *this;
  }

  friend bool operator==<>(const Point<D>& lhs, const Point<D>& rhs);

  void Print() const {
    for (size_t i = 0; i < D; i++) {
      std::cout << data_[i] << " ";
    }
    std::cout << std::endl;
  }

 private:
  std::array<double, D> data_;
};

template <size_t D>
bool operator==(const Point<D>& lhs, const Point<D>& rhs) {
  for (size_t i = 0; i < D; i++) {
    if (std::abs(lhs.data_[i] - rhs.data_[i]) > 1e-10) {
      return false;
    }
  }
  return true;
}
