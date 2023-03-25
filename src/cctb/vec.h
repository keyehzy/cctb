#pragma once

#include <vector>
#include <iostream>

template <typename T>
struct Vec {
  std::vector<T> data;
  
  Vec(T x) : data({x}) {}

  Vec(T x, T y) : data({x, y}) {}

  Vec(T x, T y, T z) : data({x, y, z}) {}

  Vec(const std::vector<T>& data) : data(data) {}

  std::size_t size() const { return data.size(); }

  T& operator[](int i) { return data[i]; }

  T operator[](int i) const { return data[i]; }

  Vec operator+(const Vec& other) const {
    Vec result(data);
    for (std::size_t i = 0; i < data.size(); i++) {
      result[i] += other[i];
    }
    return result;
  }

  Vec operator-(const Vec& other) const {
    Vec result(data);
    for (std::size_t i = 0; i < data.size(); i++) {
      result[i] -= other[i];
    }
    return result;
  }

  Vec operator*(T scalar) const {
    Vec result(data);
    for (std::size_t i = 0; i < data.size(); i++) {
      result[i] *= scalar;
    }
    return result;
  }

  bool operator==(const Vec& other) const {
    if (data.size() != other.data.size()) {
      return false;
    }
    for (std::size_t i = 0; i < data.size(); i++) {
      if (data[i] != other.data[i]) {
        return false;
      }
    }
    return true;
  }

  void Print() const {
    for (int i = 0; i < data.size(); i++) {
      std::cout << data[i] << " ";
    }
    std::cout << std::endl;
  }
};
