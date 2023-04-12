#include "Geometry/Region/Parallelogram.h"

// https://math.stackexchange.com/questions/2642918/decide-if-a-point-is-inside-parallelogram
bool Parallelogram::contains(const Point<2>& point) const {
  double u1 = u_[0];
  double u2 = u_[1];
  double v1 = v_[0];
  double v2 = v_[1];
  double det = u1 * v2 - u2 * v1;
  double p1 = point[0];
  double p2 = point[1];
  double mu = (p1 * v2 - p2 * v1) / det;
  double lambda = (p2 * u1 - p1 * u2) / det;
  return mu >= 0 && mu <= 1 && lambda >= 0 && lambda <= 1;
}

Mesh Parallelogram::MakeMesh(int n) const {
  Mesh mesh;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double x = (u_[0] * i + v_[0] * j) / static_cast<double>(n - 1);
      double y = (u_[1] * i + v_[1] * j) / static_cast<double>(n - 1);
      mesh.Add(Point<2>{x, y});
    }
  }
  return mesh;
}
