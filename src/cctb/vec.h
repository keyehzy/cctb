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

template <typename T>
struct Vec {
  std::vector<T> data;

  Vec(){};

  Vec(T x) : data({x}) {}

  Vec(T x, T y) : data({x, y}) {}

  Vec(T x, T y, T z) : data({x, y, z}) {}

  Vec(const std::vector<T>& data) : data(data) {}

  int size() const { return data.size(); }

  T& operator[](int i) { return data[i]; }

  T operator[](int i) const { return data[i]; }

  Vec operator+(const Vec& other) const {
    Vec result(data);
    for (int i = 0; i < data.size(); i++) {
      result[i] += other[i];
    }
    return result;
  }

  Vec operator-(const Vec& other) const {
    Vec result(data);
    for (int i = 0; i < data.size(); i++) {
      result[i] -= other[i];
    }
    return result;
  }

  Vec operator*(T scalar) const {
    Vec result(data);
    for (int i = 0; i < data.size(); i++) {
      result[i] *= scalar;
    }
    return result;
  }

  Vec operator/(T scalar) const {
    Vec result(data);
    for (int i = 0; i < data.size(); i++) {
      result[i] /= scalar;
    }
    return result;
  }

  bool operator==(const Vec& other) const {
    if (data.size() != other.data.size()) {
      return false;
    }
    for (int i = 0; i < data.size(); i++) {
      if (data[i] != other.data[i]) {
        return false;
      }
    }
    return true;
  }

  double dot(const Vec& other) const {
    double result = 0.0;
    for (int i = 0; i < data.size(); i++) {
      result += data[i] * other.data[i];
    }
    return result;
  }

  void Print() const {
    for (int i = 0; i < data.size(); i++) {
      std::cout << data[i] << " ";
    }
    std::cout << std::endl;
  }
};
