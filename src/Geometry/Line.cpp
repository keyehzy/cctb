#include "Geometry/Line.h"

Vec<float> Line::intercect(const Line& other) const {
  float x1 = p1_[0];
  float y1 = p1_[1];
  float x2 = p2_[0];
  float y2 = p2_[1];
  float x3 = other.p1()[0];
  float y3 = other.p1()[1];
  float x4 = other.p2()[0];
  float y4 = other.p2()[1];
  float x =
      ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) /
      ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
  float y =
      ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) /
      ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
  return Vec<float>{x, y};
}

Line Line::perpendicular_bisector() const {
  float x1 = p1_[0];
  float y1 = p1_[1];
  float x2 = p2_[0];
  float y2 = p2_[1];
  float x = (x1 + x2) / 2;
  float y = (y1 + y2) / 2;
  return Line(Vec<float>{x, y}, Vec<float>{x + -y, y + x});
}

bool Line::orthogonal(const Line& other) const {
  Vec<float> p1 = p2_ - p1_;
  Vec<float> p2 = other.p2() - other.p1();
  return std::abs(p1.dot(p2)) < 1e-6;
}

Vec<float> Line::midpoint() const {
  float x1 = p1_[0];
  float y1 = p1_[1];
  float x2 = p2_[0];
  float y2 = p2_[1];
  float x = (x1 + x2) / 2;
  float y = (y1 + y2) / 2;
  return Vec<float>{x, y};
}

void Line::Print() const {
  std::cout << "Line: (" << p1_[0] << ", " << p1_[1] << ") -> (" << p2_[0]
            << ", " << p2_[1] << ")" << std::endl;
}