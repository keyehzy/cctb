#pragma once

#include <array>
#include <iostream>
#include <vector>

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

  double& operator[](int i) { return data_[i]; }
  double operator[](int i) const { return data_[i]; }

  Point& operator=(const Point& other) {
    data_ = other.data_;
    return *this;
  }

  Point operator+(const Point& other) const {
    Point result(data_);
    for (int i = 0; i < D; i++) {
      result[i] += other[i];
    }
    return result;
  }

  Point operator-(const Point& other) const {
    Point result(data_);
    for (int i = 0; i < D; i++) {
      result[i] -= other[i];
    }
    return result;
  }

  Point operator*(double scalar) const {
    Point result(data_);
    for (int i = 0; i < D; i++) {
      result[i] *= scalar;
    }
    return result;
  }

  Point operator/(double scalar) const {
    Point result(data_);
    for (int i = 0; i < D; i++) {
      result[i] /= scalar;
    }
    return result;
  }

  bool operator==(const Point& other) const {
    for (int i = 0; i < D; i++) {
      if (data_[i] != other.data_[i]) {
        return false;
      }
    }
    return true;
  }

  double dot(const Point& other) const {
    double result = 0.0;
    for (int i = 0; i < D; i++) {
      result += data_[i] * other.data_[i];
    }
    return result;
  }

  void Print() const {
    for (int i = 0; i < D; i++) {
      std::cout << data_[i] << " ";
    }
    std::cout << std::endl;
  }

 private:
  std::array<double, D> data_;
};
