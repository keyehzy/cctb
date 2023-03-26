#include <cmath>
#include <utility>

#include "cctb/lattice.h"
#include "cctb/matrix.h"
#include "cctb/vec.h"

class LinearChain : public OneDimensionalLattice {
 public:
  LinearChain(int size) : OneDimensionalLattice(Vec<float>{1.0f * size}) {
    for (int i = 0; i < size; i++) {
      AddSite(Site(Vec<float>{(float)i}));
    }
    for (int i = 0; i < size - 1; i++) {
      AddEdge(Edge({0}, i, i + 1));
    }
    // if periodic, add the last edge
    AddEdge(Edge({1}, size - 1, 0));
  }
};

int main(void) {
  LinearChain lattice(2);
  lattice.Plot();
  lattice.AdjMatrix().Print();
  lattice.HoppingMatrix(0.5).Print();
  return 0;
}
