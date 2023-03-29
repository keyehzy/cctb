#pragma once

#include <array>
#include <iostream>
#include <vector>

template <std::size_t D>
class NewVec {
 public:
  NewVec() : data_({}){};
  NewVec(double x) : data_({x}){};
  NewVec(double x, double y) : data_({x, y}){};
  NewVec(double x, double y, double z) : data_({x, y, z}){};
  NewVec(const std::array<double, D>& a) : data_(a){};
  NewVec(const NewVec& other) : data_(other.data_){};

  int size() const { return D; }

  double& operator[](int i) { return data_[i]; }
  double operator[](int i) const { return data_[i]; }

  NewVec& operator=(const NewVec& other) {
    data_ = other.data_;
    return *this;
  }

  NewVec operator+(const NewVec& other) const {
    NewVec result(data_);
    for (int i = 0; i < D; i++) {
      result[i] += other[i];
    }
    return result;
  }

  NewVec operator-(const NewVec& other) const {
    NewVec result(data_);
    for (int i = 0; i < D; i++) {
      result[i] -= other[i];
    }
    return result;
  }

  NewVec operator*(double scalar) const {
    NewVec result(data_);
    for (int i = 0; i < D; i++) {
      result[i] *= scalar;
    }
    return result;
  }

  NewVec operator/(double scalar) const {
    NewVec result(data_);
    for (int i = 0; i < D; i++) {
      result[i] /= scalar;
    }
    return result;
  }

  bool operator==(const NewVec& other) const {
    for (int i = 0; i < D; i++) {
      if (data_[i] != other.data_[i]) {
        return false;
      }
    }
    return true;
  }

  double dot(const NewVec& other) const {
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
