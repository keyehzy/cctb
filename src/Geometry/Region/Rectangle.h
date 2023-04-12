#pragma once

#include "Geometry/Region/Parallelogram.h"

class Rectangle : public Parallelogram {
 public:
  Rectangle(const Point<2>& origin, double width, double height)
      : Parallelogram(origin, Vector<2>(width, 0), Vector<2>(0, height)) {}

  ~Rectangle() {}

  double width() const { return width_; }

  double height() const { return height_; }

 private:
  double width_;
  double height_;
};
