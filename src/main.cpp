#include <cmath>

#include "cctb/lattice.h"
#include "cctb/matrix.h"
#include "cctb/vec.h"

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

int main(void) {
  GrapheneLatticeExtended lattice;
  lattice.Plot();
  lattice.AdjMatrix().Print();
  return 0;
}
