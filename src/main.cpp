#include <cmath>
#include <fstream>
#include <utility>

#include "Geometry/Point.h"
#include "lattice.h"

class GrapheneLatticeExtendedTest : public TwoDimensionalLattice {
 public:
  GrapheneLatticeExtendedTest()
      : TwoDimensionalLattice(Vector<2>(3.0, 0), Vector<2>(0, sqrt(3.0))) {
    add_site(Point<2>(0, 0));
    add_site(Point<2>(0.5, 0.5 * sqrt(3.0)));
    add_site(Point<2>(1.5, 0.5 * sqrt(3.0)));
    add_site(Point<2>(2.0, 0));

    // intra unit cell
    add_edge(0, 1, {0, 0}, 1.0);
    add_edge(1, 2, {0, 0}, 1.0);
    add_edge(2, 3, {0, 0}, 1.0);

    // inter unit cell
    add_edge(3, 0, {1, 0}, 1.0);
    add_edge(1, 0, {0, 1}, 1.0);
    add_edge(2, 3, {0, 1}, 1.0);
  };
};

class TriangularLatticeTest : public TwoDimensionalLattice {
 public:
  TriangularLatticeTest(double a = 1.0)
      : TwoDimensionalLattice(Vector<2>(a, 0), Vector<2>(0.5 * a, 0.5 * a * sqrt(3.0))) {
    add_site(Point<2>(0, 0));
    add_edge(0, 0, {1, 0}, 1.0);

    // intra unit cell
    add_edge(0, 0, {0, 1}, 1.0);

    // inter unit cell
    add_edge(0, 0, {1, -1}, 1.0);
  }
};

class GrapheneLattice : public TwoDimensionalLattice {
 public:
  GrapheneLattice()
      : TwoDimensionalLattice(Vector<2>(1.5, 0.5 * sqrt(3.0)), Vector<2>(1.5, -0.5 * sqrt(3.0))) {
    add_site(Point<2>{0, 0});
    add_site(Point<2>{0.5, 0.5 * sqrt(3.0)});

    // intra unit cell
    add_edge(0, 1, {0, 0}, 1.0);

    // inter unit cell
    add_edge(1, 0, {1, 0}, 1.0);
    add_edge(1, 0, {1, -1}, 1.0);
  }
};

int main(void) {
  GrapheneLattice lattice;
  std::ofstream lattice_file("lattice.asy");
  lattice.Plot(PainterBackend::kAsymptote, lattice_file);
  std::ofstream bz_file("bz.asy");
  lattice.PlotBrillouinZone(PainterBackend::kAsymptote, bz_file);
  return 0;
}
