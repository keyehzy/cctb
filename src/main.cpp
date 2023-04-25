#include <cmath>
#include <complex>
#include <fstream>
#include <utility>

#include "Lattice/Lattice.h"

using namespace std::complex_literals;

class LinearChainTest : public SingleOrbitalOneDimensionalLattice {
 public:
  using NumberType = double;
  using Complex = std::complex<NumberType>;
  using Vector = Eigen::Vector<NumberType, 2>;
  using Matrix = Eigen::Matrix<Complex, 1, 1>;

  LinearChainTest(size_t size)
      : OneDimensionalLattice(VectorType(1.0 * static_cast<double>(size))) {
    Matrix onsite = Matrix::Zero();
    Matrix hopping = Matrix::Identity();

    for (size_t i = 0; i < size; i++) {
      add_site(0, VectorType(static_cast<double>(i)), onsite);
    }
    for (size_t i = 0; i < size - 1; i++) {
      add_edge(i, i + 1, {0}, hopping);
    }
    // if periodic, add the last edge
    add_edge(size - 1, 0, {1}, hopping);
  }
};

class GrapheneLatticeExtendedTest : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using ComplexType = std::complex<NumberType>;
  using VectorType = Eigen::Vector<NumberType, 2>;
  using MatrixType = Eigen::Matrix<NumberType, 1, 1>;

  GrapheneLatticeExtendedTest()
      : TwoDimensionalLattice(VectorType(3.0, 0), VectorType(0, sqrt(3.0))) {
    MatrixType onsite = MatrixType::Zero();
    MatrixType hopping = MatrixType::Identity();
    add_site(0, VectorType(0, 0), onsite);
    add_site(1, VectorType(0.5, 0.5 * sqrt(3.0)), onsite);
    add_site(2, VectorType(1.5, 0.5 * sqrt(3.0)), onsite);
    add_site(3, VectorType(2.0, 0), onsite);

    add_edge(0, 1, {0, 0}, hopping);
    add_edge(1, 2, {0, 0}, hopping);
    add_edge(2, 3, {0, 0}, hopping);

    add_edge(3, 0, {1, 0}, hopping);
    add_edge(1, 0, {0, 1}, hopping);
    add_edge(2, 3, {0, 1}, hopping);
  }
};

class TriangularLatticeTest : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using ComplexType = std::complex<NumberType>;
  using VectorType = Eigen::Vector<NumberType, 2>;
  using MatrixType = Eigen::Matrix<NumberType, 1, 1>;

  TriangularLatticeTest()
      : TwoDimensionalLattice(VectorType(1, 0), VectorType(0.5, 0.5 * sqrt(3.0))) {
    MatrixType onsite = MatrixType::Zero();
    MatrixType hopping = MatrixType::Identity();
    add_site(0, VectorType(0, 0), onsite);
    add_edge(0, 0, {1, 0}, hopping);
    add_edge(0, 0, {0, 1}, hopping);
    add_edge(0, 0, {1, -1}, hopping);
  }
};

class GrapheneLatticeTest : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using ComplexType = std::complex<NumberType>;
  using VectorType = Eigen::Vector<NumberType, 2>;
  using MatrixType = Eigen::Matrix<NumberType, 1, 1>;

  GrapheneLatticeTest()
      : TwoDimensionalLattice(VectorType(1.5, 0.5 * sqrt(3.0)), VectorType(1.5, -0.5 * sqrt(3.0))) {
    MatrixType onsite = MatrixType::Zero();
    MatrixType hopping = MatrixType::Identity();
    add_site(0, VectorType(0, 0), onsite);
    add_site(1, VectorType(0.5, 0.5 * sqrt(3.0)), onsite);
    add_edge(0, 1, {0, 0}, hopping);
    add_edge(1, 0, {1, 0}, hopping);
    add_edge(1, 0, {1, -1}, hopping);
  }
};

class SquareLatticeExtended : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using ComplexType = std::complex<NumberType>;
  using VectorType = Eigen::Vector<NumberType, 2>;
  using MatrixType = Eigen::Matrix<NumberType, 1, 1>;

  SquareLatticeExtended() : TwoDimensionalLattice(VectorType(2.0, 0), VectorType(0, 2.0)) {
    MatrixType onsite = MatrixType::Zero();
    MatrixType hopping = MatrixType::Identity();

    add_site(0, VectorType(0, 0), onsite);
    add_site(1, VectorType(1, 0), onsite);
    add_site(2, VectorType(0, 1), onsite);
    add_site(3, VectorType(1, 1), onsite);

    add_edge(0, 1, {0, 0}, hopping);
    add_edge(0, 2, {0, 0}, hopping);
    add_edge(1, 3, {0, 0}, hopping);
    add_edge(2, 3, {0, 0}, hopping);

    add_edge(1, 0, {1, 0}, hopping);
    add_edge(2, 0, {0, 1}, hopping);
    add_edge(3, 1, {0, 1}, hopping);
    add_edge(3, 2, {1, 0}, hopping);
  }
};

class LiebLattice : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using Complex = std::complex<NumberType>;
  using Vector = Eigen::Vector<NumberType, 2>;
  using Matrix = Eigen::Matrix<Complex, 1, 1>;

  LiebLattice() : TwoDimensionalLattice(Vector(2.0, 0.0), Vector(0.0, 2.0)) {
    Matrix onsite = Matrix::Zero();
    Matrix hopping = Matrix::Identity();

    add_site(0, Vector(0.0, 0.0), onsite);
    add_site(1, Vector(1.0, 0.0), onsite);
    add_site(2, Vector(0.0, 1.0), onsite);

    add_edge(0, 1, {0, 0}, hopping);
    add_edge(0, 2, {0, 0}, hopping);

    add_edge(1, 0, {1, 0}, hopping);
    add_edge(2, 0, {0, 1}, hopping);
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
  SquareLatticeExtended lattice;
  // lattice.statistics();
  // std::cout << lattice.hopping_matrix({0.5, 0.8}) << std::endl;

  std::ofstream lattice_file("lattice.tex");
  lattice.Plot(PainterBackend::kTikz, lattice_file);
  // std::ofstream bz_file("bz.asy");
  // lattice.PlotBrillouinZone(PainterBackend::kAsymptote, bz_file);
  std::ofstream band_file("band_new.tex");
  lattice.PlotBandStructure(band_file);
  return 0;
}
