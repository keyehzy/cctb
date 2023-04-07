#pragma once

#include "Geometry/Region/Region.h"

class Rectangle : public Region<2> {
 public:
  Rectangle(){};
  Rectangle(const Point<2>& origin, double width, double height)
      : origin_(origin), width_(width), height_(height){};
  virtual ~Rectangle(){};

  Point<2> origin() const override { return origin_; }
  bool contains(const Point<2>& point) const override;
  std::vector<Point<2>> grid(size_t n) const override;
  double width() const { return width_; }
  double height() const { return height_; }

 private:
  Point<2> origin_;
  double width_;
  double height_;
};
