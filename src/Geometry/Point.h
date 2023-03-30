#pragma once

#include <array>
#include <iostream>

#include "Geometry/Vector.h"

template <std::size_t D>
class Vector;
template <std::size_t D>
class Point;

template <std::size_t D>
bool operator==(const Point<D>& lhs, const Point<D>& rhs);

template <std::size_t D>
class Point {
 public:
  Point() : data_({}){};
  Point(double x) : data_({x}){};
  Point(double x, double y) : data_({x, y}){};
  Point(double x, double y, double z) : data_({x, y, z}){};
  Point(const std::array<double, D>& a) : data_(a){};
  Point(const Point& other) : data_(other.data_){};

  int size() const { return D; }

  Point translated(const Vector<D>& v) const {
    Point result;
    for (int i = 0; i < D; i++) {
      result.data_[i] = data_[i] + v[i];
    }
    return result;
  }

  double& operator[](int i) { return data_[i]; }
  double operator[](int i) const { return data_[i]; }

  Point& operator=(const Point& other) {
    data_ = other.data_;
    return *this;
  }

  friend bool operator==<>(const Point<D>& lhs, const Point<D>& rhs);

  void Print() const {
    for (int i = 0; i < D; i++) {
      std::cout << data_[i] << " ";
    }
    std::cout << std::endl;
  }

 private:
  std::array<double, D> data_;
};

template <std::size_t D>
bool operator==(const Point<D>& lhs, const Point<D>& rhs) {
  for (int i = 0; i < D; i++) {
    if (lhs.data_[i] != rhs.data_[i]) {
      return false;
    }
  }
  return true;
}
