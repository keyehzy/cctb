#pragma once

#include <array>
#include <cmath>
#include <iostream>

#include "Geometry/Point.h"

template <std::size_t D>
class Vector;
template <std::size_t D>
class Point;

template <std::size_t D>
bool operator==(const Vector<D>& lhs, const Vector<D>& rhs);
template <std::size_t D>
Vector<D> operator+(const Vector<D>& lhs, const Vector<D>& rhs);
template <std::size_t D>
Vector<D> operator-(const Vector<D>& lhs, const Vector<D>& rhs);
template <std::size_t D>
Vector<D> operator*(const Vector<D>& lhs, double rhs);
template <std::size_t D>
Vector<D> operator*(double lhs, const Vector<D>& rhs);
template <std::size_t D>
Vector<D> operator/(const Vector<D>& lhs, double rhs);
template <std::size_t D>
Vector<D> operator/(double lhs, const Vector<D>& rhs);

template <std::size_t D>
class Vector {
 public:
  Vector() : data_({}){};
  Vector(double x) : data_({x}){};
  Vector(double x, double y) : data_({x, y}){};
  Vector(double x, double y, double z) : data_({x, y, z}){};
  Vector(const std::array<double, D>& a) : data_(a){};
  Vector(const Vector& other) : data_(other.data_){};

  Vector(const Point<D>& p1, const Point<D>& p2) {
    for (int i = 0; i < D; i++) {
      data_[i] = p2[i] - p1[i];
    }
  }

  int size() const { return D; }

  double& operator[](int i) { return data_[i]; }
  double operator[](int i) const { return data_[i]; }

  double dot(const Vector& other) const {
    double result = 0;
    for (int i = 0; i < D; i++) {
      result += data_[i] * other.data_[i];
    }
    return result;
  }

  double norm() const {
    double result = 0;
    for (int i = 0; i < D; i++) {
      result += data_[i] * data_[i];
    }
    return std::sqrt(result);
  }

  void scale(double alpha) {
    for (int i = 0; i < D; i++) {
      data_[i] *= alpha;
    }
  }

  void normalize() {
    double n = norm();
    for (int i = 0; i < D; i++) {
      data_[i] /= n;
    }
  }

  Vector& operator=(const Vector& other) {
    data_ = other.data_;
    return *this;
  }

  Vector& operator+=(const Vector& other) {
    for (int i = 0; i < D; i++) {
      data_[i] += other.data_[i];
    }
    return *this;
  }

  Vector& operator-=(const Vector& other) {
    for (int i = 0; i < D; i++) {
      data_[i] -= other.data_[i];
    }
    return *this;
  }

  Vector& operator*=(double alpha) {
    for (int i = 0; i < D; i++) {
      data_[i] *= alpha;
    }
    return *this;
  }

  Vector& operator/=(double alpha) {
    for (int i = 0; i < D; i++) {
      data_[i] /= alpha;
    }
    return *this;
  }

  friend bool operator==<>(const Vector& lhs, const Vector& rhs);
  friend Vector operator+<>(const Vector& lhs, const Vector& rhs);
  friend Vector operator-<>(const Vector& lhs, const Vector& rhs);
  friend Vector operator*<>(const Vector& lhs, double rhs);
  friend Vector operator*<>(double lhs, const Vector& rhs);
  friend Vector operator/<>(const Vector& lhs, double rhs);
  friend Vector operator/<>(double lhs, const Vector& rhs);

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
bool operator==(const Vector<D>& lhs, const Vector<D>& rhs) {
  for (int i = 0; i < D; i++) {
    if (lhs[i] != rhs[i]) {
      return false;
    }
  }
  return true;
}

template <std::size_t D>
Vector<D> operator+(const Vector<D>& lhs, const Vector<D>& rhs) {
  Vector<D> result;
  for (int i = 0; i < D; i++) {
    result[i] = lhs[i] + rhs[i];
  }
  return result;
}

template <std::size_t D>
Vector<D> operator-(const Vector<D>& lhs, const Vector<D>& rhs) {
  Vector<D> result;
  for (int i = 0; i < D; i++) {
    result[i] = lhs[i] - rhs[i];
  }
  return result;
}

template <std::size_t D>
Vector<D> operator*(const Vector<D>& lhs, double rhs) {
  Vector<D> result;
  for (int i = 0; i < D; i++) {
    result[i] = lhs[i] * rhs;
  }
  return result;
}

template <std::size_t D>
Vector<D> operator*(double lhs, const Vector<D>& rhs) {
  return rhs * lhs;
}

template <std::size_t D>
Vector<D> operator/(const Vector<D>& lhs, double rhs) {
  return lhs * (1.0 / rhs);
}

template <std::size_t D>
Vector<D> operator/(double lhs, const Vector<D>& rhs) {
  return rhs / lhs;
}
