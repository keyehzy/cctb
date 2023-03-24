#include "cctb/lattice.h"

int main(void)
{
    OneDimensionalLattice lattice(Vec2f(0, 1), Vec2f(1, 0));
    lattice.AddSite(Vec2f(0, 0));
    lattice.AddEdge(Edge({1}, 0, 0));
    return 0;
}