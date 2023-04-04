#pragma once

#include "LinearAlgebra/NBuffer.h"

template <typename Lhs, typename Rhs, typename Op>
struct LazyOp {
  Lhs const& lhs_;
  Rhs const& rhs_;
  Op const& op_;

  LazyOp(Lhs const& lhs, Rhs const& rhs, Op const& op) : lhs_(lhs), rhs_(rhs), op_(op) {}

  template <typename T>
  operator T() const {
    return op_(lhs_, rhs_);
  }
};

template <typename Lhs, typename Rhs>
struct AddOp : public LazyOp<Lhs, Rhs, AddOp<Lhs, Rhs>> {
  AddOp(Lhs const& lhs, Rhs const& rhs) : LazyOp<Lhs, Rhs, AddOp<Lhs, Rhs>>(lhs, rhs, *this) {}

  template <typename T>
  T operator()(T const& lhs, T const& rhs) const {
    return lhs + rhs;
  }
};

template <typename Lhs, typename Rhs>
struct SubOp : public LazyOp<Lhs, Rhs, SubOp<Lhs, Rhs>> {
  SubOp(Lhs const& lhs, Rhs const& rhs) : LazyOp<Lhs, Rhs, SubOp<Lhs, Rhs>>(lhs, rhs, *this) {}

  template <typename T>
  T operator()(T const& lhs, T const& rhs) const {
    return lhs - rhs;
  }
};

template <typename Lhs, typename Rhs>
struct MulOp : public LazyOp<Lhs, Rhs, MulOp<Lhs, Rhs>> {
  MulOp(Lhs const& lhs, Rhs const& rhs) : LazyOp<Lhs, Rhs, MulOp<Lhs, Rhs>>(lhs, rhs, *this) {}

  template <typename T>
  T operator()(T const& lhs, T const& rhs) const {
    return lhs * rhs;
  }
};

template <typename Lhs, typename Rhs>
struct DivOp : public LazyOp<Lhs, Rhs, DivOp<Lhs, Rhs>> {
  DivOp(Lhs const& lhs, Rhs const& rhs) : LazyOp<Lhs, Rhs, DivOp<Lhs, Rhs>>(lhs, rhs, *this) {}

  template <typename T>
  T operator()(T const& lhs, T const& rhs) const {
    return lhs / rhs;
  }
};

template <typename Lhs, typename Rhs, typename Rhs2>
struct FmaOp : public LazyOp<Lhs, Rhs, FmaOp<Lhs, Rhs, Rhs2>> {
  Rhs2 const& rhs2_;

  FmaOp(Lhs const& lhs, Rhs const& rhs, Rhs2 const& rhs2)
      : LazyOp<Lhs, Rhs, FmaOp<Lhs, Rhs, Rhs2>>(lhs, rhs, *this), rhs2_(rhs2) {}

  template <typename T>
  T operator()(T const& lhs, T const& rhs) const {
    return lhs * rhs + rhs2_;
  }
};

template <typename Lhs, typename Rhs>
auto operator+(Lhs const& lhs, Rhs const& rhs) {
  return AddOp<Lhs, Rhs>(lhs, rhs);
}

template <typename Lhs, typename Rhs>
auto operator-(Lhs const& lhs, Rhs const& rhs) {
  return SubOp<Lhs, Rhs>(lhs, rhs);
}

template <typename Lhs, typename Rhs>
auto operator*(Lhs const& lhs, Rhs const& rhs) {
  return MulOp<Lhs, Rhs>(lhs, rhs);
}

template <typename Lhs, typename Rhs>
auto operator/(Lhs const& lhs, Rhs const& rhs) {
  return DivOp<Lhs, Rhs>(lhs, rhs);
}

template <typename Lhs, typename Rhs, typename Rhs2>
auto fma(Lhs const& lhs, Rhs const& rhs, Rhs2 const& rhs2) {
  return FmaOp<Lhs, Rhs, Rhs2>(lhs, rhs, rhs2);
}
