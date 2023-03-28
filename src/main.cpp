#include <cmath>
#include <fstream>
#include <utility>

#include "cctb/lattice.h"
#include "cctb/matrix.h"
#include "cctb/vec.h"

class GrapheneLatticeExtendedTest : public TwoDimensionalLattice {
 public:
  GrapheneLatticeExtendedTest()
      : TwoDimensionalLattice(Vec<double>(3.0, 0), Vec<double>(0, sqrt(3.0))) {
    AddSite(Site(Vec<double>{0, 0}));
    AddSite(Site(Vec<double>{0.5, 0.5 * sqrt(3.0)}));
    AddSite(Site(Vec<double>{1.5, 0.5 * sqrt(3.0)}));
    AddSite(Site(Vec<double>{2.0, 0}));
    AddEdge(Edge({0, 0}, 0, 1));
    AddEdge(Edge({0, 0}, 1, 2));
    AddEdge(Edge({0, 0}, 2, 3));

    AddEdge(Edge({1, 0}, 3, 0));
    AddEdge(Edge({0, 1}, 1, 0));
    AddEdge(Edge({0, 1}, 2, 3));
  };
};

class TriangularLatticeTest : public TwoDimensionalLattice {
 public:
  TriangularLatticeTest(double a = 1.0)
      : TwoDimensionalLattice(Vec<double>(a, 0), Vec<double>(0.5 * a, 0.5 * a * sqrt(3.0))) {
    AddSite(Site(Vec<double>{0, 0}));
    AddEdge(Edge({1, 0}, 0, 0));
    AddEdge(Edge({0, 1}, 0, 0));
    AddEdge(Edge({1, -1}, 0, 0));
  }
};

class GrapheneLattice : public TwoDimensionalLattice {
 public:
  GrapheneLattice()
      : TwoDimensionalLattice(Vec<double>(1.5, 0.5 * sqrt(3.0)),
                              Vec<double>(1.5, -0.5 * sqrt(3.0))) {
    AddSite(Site(Vec<double>{0, 0}));
    AddSite(Site(Vec<double>{0.5, 0.5 * sqrt(3.0)}));
    AddEdge(Edge({0, 0}, 0, 1));
    AddEdge(Edge({1, 0}, 1, 0));
    AddEdge(Edge({1, -1}, 1, 0));
  }
};

int main(void) {
  TriangularLatticeTest lattice;
  std::ofstream lattice_file("lattice.tex");
  lattice.Plot(PainterBackend::kTikz, lattice_file);
  // lattice.AdjMatrix().Print();
  // lattice.HoppingMatrix(Vec<double>{0.5, 0.8}).Print();

  std::ofstream bz_file("bz.tex");
  lattice.PlotBrillouinZone(PainterBackend::kTikz, bz_file);
  return 0;
}
