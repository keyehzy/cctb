#include <cmath>
#include <complex>
#include <fstream>
#include <utility>

#include "Lattice/Lattice.h"

using namespace std::complex_literals;

/*
  class LinearChainTest : public OneDimensionalLattice {
  public:
  LinearChainTest(size_t size) : OneDimensionalLattice(Vector<1>(1.0 * static_cast<double>(size))) {
  for (size_t i = 0; i < size; i++) {
  add_site(Point<1>(static_cast<double>(i)));
  }
  for (size_t i = 0; i < size - 1; i++) {
  add_edge(i, i + 1, {0}, 1.0);
  }
  // if periodic, add the last edge
  add_edge(size - 1, 0, {1}, 1.0);
  }
  };

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
  }
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
*/
class GrapheneLattice : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using Complex = std::complex<NumberType>;
  using Vector = Eigen::Vector<NumberType, 2>;
  using Matrix = Eigen::Matrix<Complex, 1, 1>;

  GrapheneLattice()
      : TwoDimensionalLattice(Vector(1.5, 0.5 * sqrt(3.0)), Vector(1.5, -0.5 * sqrt(3.0))) {
    Matrix onsite = Matrix::Zero();
    Matrix hopping = Matrix::Identity();

    add_site(0, Vector(0, 0), onsite);
    add_site(0, Vector(0.5, 0.5 * sqrt(3)), onsite);

    // intra unit cell
    add_edge(0, 1, {0, 0}, hopping);

    // inter unit cell
    add_edge(1, 0, {1, 0}, hopping);
    add_edge(1, 0, {1, -1}, hopping);
  }
};

class GrapheneLatticeRashba : public TwoDimensionalLattice<2> {
 public:
  using NumberType = double;
  using Complex = std::complex<NumberType>;
  using Vector = Eigen::Vector<NumberType, 2>;
  using Matrix = Eigen::Matrix<Complex, 2, 2>;

  GrapheneLatticeRashba()
      : TwoDimensionalLattice(Vector(1.5, 0.5 * sqrt(3.0)), Vector(1.5, -0.5 * sqrt(3.0))) {
    Matrix onsite = Matrix::Zero();
    Matrix hopping = Matrix::Identity();
    double rashba = 0.1;

    Matrix sigma_0{{1, 0}, {0, 1}};
    Matrix sigma_x{{0, 1}, {1, 0}};
    Matrix sigma_y{{0, -1i}, {1i, 0}};
    Matrix sigma_z{{1, 0}, {0, -1}};

    Vector d1{0.5, 0.5 * sqrt(3)};
    Vector d2{1.5, -0.5 * sqrt(3)};
    Vector d3{-1.0, 0.0};

    Matrix t1 = hopping + rashba * (sigma_x * d1[1] - sigma_y * d1[0]);
    Matrix t2 = hopping + rashba * (sigma_x * d2[1] - sigma_y * d2[0]);
    Matrix t3 = hopping + rashba * (sigma_x * d3[1] - sigma_y * d3[0]);

    add_site(0, Vector(0, 0), onsite);
    add_site(0, Vector(0.5, 0.5 * sqrt(3)), onsite);

    // intra unit cell
    add_edge(0, 1, {0, 0}, t1);

    // inter unit cell
    add_edge(1, 0, {1, 0}, t3);
    add_edge(1, 0, {1, -1}, t2);
  }
};

class GrapheneLatticeExtended : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using Complex = std::complex<NumberType>;
  using Vector = Eigen::Vector<NumberType, 2>;
  using Matrix = Eigen::Matrix<Complex, 1, 1>;

  GrapheneLatticeExtended() : TwoDimensionalLattice(Vector(3.0, 0), Vector(0, sqrt(3.0))) {
    Matrix onsite = Matrix::Zero();
    Matrix hopping = Matrix::Identity();

    add_site(0, Vector(0, 0), onsite);
    add_site(1, Vector(0.5, 0.5 * sqrt(3.0)), onsite);
    add_site(2, Vector(1.5, 0.5 * sqrt(3.0)), onsite);
    add_site(3, Vector(2.0, 0), onsite);

    // intra unit cell
    add_edge(0, 1, {0, 0}, hopping);
    add_edge(1, 2, {0, 0}, hopping);
    add_edge(2, 3, {0, 0}, hopping);

    // inter unit cell
    add_edge(3, 0, {1, 0}, hopping);
    add_edge(1, 0, {0, 1}, hopping);
    add_edge(2, 3, {0, 1}, hopping);
  }
};

int main(void) {
  GrapheneLatticeRashba lattice;
  lattice.statistics();
  // std::cout << lattice.hopping_matrix({0.5, 0.8}) << std::endl;

  // std::ofstream lattice_file("lattice.asy");
  // lattice.Plot(PainterBackend::kAsymptote, lattice_file);
  // std::ofstream bz_file("bz.asy");
  // lattice.PlotBrillouinZone(PainterBackend::kAsymptote, bz_file);
  // std::ofstream band_file("band_new.tex");
  // lattice.PlotBandStructure(band_file);
  return 0;
}
