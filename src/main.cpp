#include <cmath>
#include <fstream>
#include <utility>

#include "cctb/lattice.h"
#include "cctb/matrix.h"
#include "cctb/vec.h"

class GrapheneLattice : public TwoDimensionalLattice {
 public:
  GrapheneLattice()
      : TwoDimensionalLattice(Vec<float>(1.5f, 0.5f * sqrtf(3.0f)),
                              Vec<float>(1.5f, -0.5f * sqrtf(3.0f))) {
    AddSite(Site(Vec<float>{0, 0}));
    AddSite(Site(Vec<float>{0.5, 0.5f * sqrtf(3.0f)}));
    AddEdge(Edge({0, 0}, 0, 1));
    AddEdge(Edge({1, 0}, 1, 0));
    AddEdge(Edge({1, -1}, 1, 0));
  }
};

int main(void) {
  GrapheneLattice lattice;
  std::ofstream lattice_file("lattice.tex");
  lattice.Plot(PainterBackend::kTikz, lattice_file);
  // lattice.AdjMatrix().Print();
  // lattice.HoppingMatrix(Vec<float>{0.5, 0.8}).Print();

  std::ofstream bz_file("bz.tex");
  lattice.PlotBrillouinZone(PainterBackend::kTikz, bz_file);
  return 0;
}
