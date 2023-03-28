#include <cmath>
#include <fstream>
#include <utility>

#include "cctb/lattice.h"
#include "cctb/matrix.h"
#include "cctb/vec.h"

class GrapheneLatticeExtendedTest : public TwoDimensionalLattice {
 public:
  GrapheneLatticeExtendedTest()
      : TwoDimensionalLattice(Vec<double>(3.0f, 0), Vec<double>(0, sqrtf(3.0f))) {
    AddSite(Site(Vec<double>{0, 0}));
    AddSite(Site(Vec<double>{0.5f, 0.5f * sqrtf(3.0f)}));
    AddSite(Site(Vec<double>{1.5f, 0.5f * sqrtf(3.0f)}));
    AddSite(Site(Vec<double>{2.0f, 0}));
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
  TriangularLatticeTest(double a = 1.0f)
      : TwoDimensionalLattice(Vec<double>(a, 0), Vec<double>(0.5f * a, 0.5f * a * sqrtf(3.0f))) {
    AddSite(Site(Vec<double>{0, 0}));
    AddEdge(Edge({1, 0}, 0, 0));
    AddEdge(Edge({0, 1}, 0, 0));
    AddEdge(Edge({1, -1}, 0, 0));
  }
};

class GrapheneLattice : public TwoDimensionalLattice {
 public:
  GrapheneLattice()
      : TwoDimensionalLattice(Vec<double>(1.5f, 0.5f * sqrtf(3.0f)),
                              Vec<double>(1.5f, -0.5f * sqrtf(3.0f))) {
    AddSite(Site(Vec<double>{0, 0}));
    AddSite(Site(Vec<double>{0.5, 0.5f * sqrtf(3.0f)}));
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
