#include <cmath>

#include "cctb/lattice.h"
#include "cctb/matrix.h"
#include "cctb/vec.h"

class SquareLattice : public TwoDimensionalLattice {
 public:
  SquareLattice(float d) : TwoDimensionalLattice(Vec<float>{d, 0}, Vec<float>{0, d}) {
    AddSite(Site(Vec<float>{0, 0}));
    AddEdge(Edge({0, 1}, 0, 0));
    AddEdge(Edge({1, 0}, 0, 0));
  }
};

class GrapheneLatticeTest : public TwoDimensionalLattice {
 public:
  GrapheneLatticeTest(float a = 0.24595f, float a_cc = 0.142f) : TwoDimensionalLattice(Vec<float>(a, 0), Vec<float>(0.5f * a, 0.5f * a * sqrtf(3.0f))) {
    AddSite(Site(Vec<float>{0, -0.5f * a_cc}));
    AddSite(Site(Vec<float>{0, 0.5f * a_cc}));
    AddEdge(Edge({0, 0}, 0, 1));
    AddEdge(Edge({1, -1}, 0, 1));
    AddEdge(Edge({0, -1}, 0, 1));
  }
};
class GrapheneLattice : public TwoDimensionalLattice {
 public:
  GrapheneLattice() : TwoDimensionalLattice(Vec<float>(1.5f, 0.5f * sqrtf(3.0f)), Vec<float>(1.5f, -0.5f * sqrtf(3.0f))) {
    AddSite(Site(Vec<float>{0, 0}));
    AddSite(Site(Vec<float>{0.5, 0.5f * sqrtf(3.0f)}));
    AddEdge(Edge({0, 0}, 0, 1));
    AddEdge(Edge({1, 0}, 1, 0));
    AddEdge(Edge({1, -1}, 1, 0));
  }
};

class GrapheneLatticeExtended : public TwoDimensionalLattice {
 public:
  GrapheneLatticeExtended() : TwoDimensionalLattice(Vec<float>(3.0f, 0), Vec<float>(0, sqrtf(3.0f))) {
    AddSite(Site(Vec<float>{0, 0}));
    AddSite(Site(Vec<float>{0.5f, 0.5f * sqrtf(3.0f)}));
    AddSite(Site(Vec<float>{1.5f, 0.5f * sqrtf(3.0f)}));
    AddSite(Site(Vec<float>{2.0f, 0}));    
    AddEdge(Edge({0, 0}, 0, 1));
    AddEdge(Edge({0, 0}, 1, 2));
    AddEdge(Edge({0, 0}, 2, 3));

    AddEdge(Edge({1, 0}, 3, 0));
    AddEdge(Edge({0, 1}, 1, 0));
    AddEdge(Edge({0, 1}, 2, 3));

  };
};

class TriangularLattice : public TwoDimensionalLattice {
 public:
  TriangularLattice(float a = 1.0f) : TwoDimensionalLattice(Vec<float>(a, 0), Vec<float>(0.5f * a, 0.5f * a * sqrtf(3.0f))) {
    AddSite(Site(Vec<float>{0, 0}));
    AddEdge(Edge({1, 0}, 0, 0));
    AddEdge(Edge({0, 1}, 0, 0));
    AddEdge(Edge({1, -1}, 0, 0));
  }
};

class KagomeLattice : public TwoDimensionalLattice {
 public:
  KagomeLattice() : TwoDimensionalLattice(Vec<float>(2, 0), Vec<float>(1.0f, sqrtf(3.0f))) {
    AddSite(Site(Vec<float>{0, 0}));
    AddSite(Site(Vec<float>{1.0f, 0}));
    AddSite(Site(Vec<float>{0.5f, 0.5f * sqrtf(3.0f)}));
    AddEdge(Edge({0, 0}, 0, 1));
    AddEdge(Edge({0, 0}, 1, 2));
    AddEdge(Edge({0, 0}, 2, 0));
    AddEdge(Edge({1, 0}, 1, 0));
    AddEdge(Edge({0, 1}, 2, 0));
    AddEdge(Edge({1, -1}, 1, 2));
  }
};

int main(void) {
  GrapheneLatticeExtended lattice;
  lattice.PlotSites();
  //lattice.AdjMatrix().print();
  return 0;
}
