#pragma once
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <complex>

#include "Geometry/Vector.h"

template <std::size_t D>
struct ApproxEqualVec : Catch::Matchers::MatcherGenericBase {
 public:
  ApproxEqualVec(Vector<D> v) : v_(v) {}

  bool match(Vector<D> const &v) const {
    if (v.size() != v_.size()) {
      return false;
    }
    for (size_t i = 0; i < v.size(); i++) {
      if (std::abs(v[i] - v_[i]) > 1e-10) {
        return false;
      }
    }
    return true;
  }

  std::string describe() const override { return "ApproxEqualVec"; }

 private:
  Vector<D> v_;
};

struct ApproxEqualComplex : Catch::Matchers::MatcherGenericBase {
 public:
  ApproxEqualComplex(std::complex<double> v) : v_(v) {}

  bool match(std::complex<double> const &v) const {
    return std::abs(v.real() - v_.real()) < 1e-10 && std::abs(v.imag() - v_.imag()) < 1e-10;
  }

  std::string describe() const override {
    return "ApproxEqualComplex:" + std::to_string(v_.real()) + "+" + std::to_string(v_.imag()) +
           "i";
  }

 private:
  std::complex<double> v_;
};
