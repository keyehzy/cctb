#pragma once

/* Float point helpers */

#include <cmath>

#define ALWAYS_INLINE __attribute__((always_inline)) inline

constexpr float float_epsilon = 1e-6f;
constexpr double double_epsilon = 1e-6;

ALWAYS_INLINE constexpr bool fp_eq(float a, float b) { return std::abs(a - b) < float_epsilon; }
ALWAYS_INLINE constexpr bool fp_eq(double a, double b) { return std::abs(a - b) < double_epsilon; }
