#pragma once

struct Vec2f {
  float x;
  float y;

  Vec2f(float x, float y) : x(x), y(y) {}

  Vec2f operator+(const Vec2f& other) {
    return Vec2f(x + other.x, y + other.y);
  }

  Vec2f operator-(const Vec2f& other) {
    return Vec2f(x - other.x, y - other.y);
  }

  Vec2f operator*(const float& other) { return Vec2f(x * other, y * other); }

  bool operator==(const Vec2f& other) { return x == other.x && y == other.y; }
};