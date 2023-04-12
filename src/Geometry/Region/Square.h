#pragma once

#include "Geometry/Region/Rectangle.h"

class Square : public Rectangle {
 public:
  Square(const Point<2>& origin, double length) : Rectangle(origin, length, length) {}

  ~Square() {}

  double length() const { return m_length; }

 private:
  double m_length;
};
