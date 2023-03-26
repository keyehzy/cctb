![Build Status](https://github.com/keyehzy/cctb/actions/workflows/cmake.yml/badge.svg)
![Clang Format](https://github.com/keyehzy/cctb/actions/workflows/clang-format.yml/badge.svg)

# CCTB --- A tightbinding software written in C++

This is a software package for simulating tightbinding models. It aims to be a simple and easy to use, making adding new models straightforward. Once you have a model, you can use it to calculate the band structure, density of states, and other properties.

## Building:

To build this project, you need to have a C++ compiler that supports C++11 and CMake. The following instructions are for Linux and macOS.

```bash
cmake -S . -B build
cmake --build build
```

The executable will be in `build/src`.

## Usage:

```bash
./build/src/cctb --help
```

## Examples:

Right now you can set up 1D and 2D lattices and graph their crystal structure. The following example shows how to set up a Graphene lattice and graph it.

```c++
#include "cctb/lattice.h"
#include "cctb/matrix.h"
#include "cctb/vec.h"

// https://arxiv.org/abs/0709.1163
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

int main(void) {
  GrapheneLattice lattice;
  std::ofstream lattice_file("lattice.tex");
  lattice.Plot(PainterBackend::kTikz, lattice_file);
  lattice.AdjMatrix().Print();
  lattice.HoppingMatrix({0.5, 0.8}).Print();
  return 0;
}
```

There are two graphing backends: `Tikz` and `Asymptote`. The each backend will generate a file that can be compiled to a PDF or any other format supported by the backend using `pdflatex` or `asy`, respectively. In the example above, we use the `Tikz` backend to generate a file called `lattice.tex`. We can then compile it to a PDF using `latexmk`:

```bash
latexmk -pdf lattice.tex
```

which will generate the following image:

![Graphene lattice](examples/graphene.pdf)


Additionally, in the code above, we also print the adjacency matrix and the hopping matrix for the lattice at a given k-point:

```bash
0 1 
1 0 
(0,0) (2.36864,-0.0986967) 
(2.36864,0.0986967) (0,0) 
```

## License:

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Related projects:

- [pybinding](https://github.com/dean0x7d/pybinding)