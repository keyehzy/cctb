#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <complex>

#include "Geometry/Point.h"
#include "Lattice/Lattice.h"
#include "LinearAlgebra/LapackImpl.h"
#include "TestMatchers.h"

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

TEST_CASE("LinearChain", "[lattice]") {
  LinearChainTest lattice(2);
  REQUIRE(lattice.size() == 2);
  REQUIRE(lattice.site(0).position == Point<1>(0));
  REQUIRE(lattice.site(0).size() == 1);
  REQUIRE(lattice.site(0).edge(0).dst == 1);
  REQUIRE(lattice.site(0).edge(0).offset == std::array<int, 1>{0});
  REQUIRE_THAT(lattice.site(0).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(1).position == Point<1>(1));
  REQUIRE(lattice.site(1).size() == 1);
  REQUIRE(lattice.site(1).edge(0).dst == 0);
  REQUIRE(lattice.site(1).edge(0).offset == std::array<int, 1>{1});
  REQUIRE_THAT(lattice.site(1).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows() == 2);
  REQUIRE(adj_matrix.cols() == 2);
  REQUIRE(adj_matrix(0, 0) == 0);
  REQUIRE(adj_matrix(0, 1) == 1);
  REQUIRE(adj_matrix(1, 0) == 1);
  REQUIRE(adj_matrix(1, 1) == 0);

  Matrix<std::complex<double>> hopping_matrix = lattice.HoppingMatrix(0.5);
  REQUIRE(hopping_matrix.rows() == 2);
  REQUIRE(hopping_matrix.cols() == 2);
  REQUIRE(hopping_matrix(0, 0) == std::complex<double>(0, 0));
  REQUIRE(hopping_matrix(0, 1) == std::complex<double>(2.0 * cos(0.5), 0));
  REQUIRE(hopping_matrix(1, 0) == std::complex<double>(2.0 * cos(0.5), 0));
  REQUIRE(hopping_matrix(1, 1) == std::complex<double>(0, 0));

  NumericArray<double> w(2);
  Matrix<std::complex<double>> v(2, 2);
  diagonalize_hermitian(hopping_matrix, w, v);
  REQUIRE_THAT(w[0], Catch::Matchers::WithinAbs(-2.0 * cos(0.5), 1e-10));
  REQUIRE_THAT(w[1], Catch::Matchers::WithinAbs(2.0 * cos(0.5), 1e-10));
}

class SquareLatticeTest : public TwoDimensionalLattice {
 public:
  SquareLatticeTest() : TwoDimensionalLattice(Vector<2>{1.0, 0}, Vector<2>(0, 1.0)) {
    add_site(Point<2>(0, 0));
    add_edge(0, 0, {0, 1}, 1.0);
    add_edge(0, 0, {1, 0}, 1.0);
  }
};

TEST_CASE("SquareLattice", "[lattice]") {
  SquareLatticeTest lattice;
  REQUIRE(lattice.size() == 1);

  REQUIRE(lattice.site(0).position == Point<2>(0, 0));
  REQUIRE(lattice.site(0).size() == 2);
  REQUIRE(lattice.site(0).edge(0).dst == 0);
  REQUIRE(lattice.site(0).edge(0).offset == std::array<int, 2>{0, 1});
  REQUIRE_THAT(lattice.site(0).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(0).edge(1).dst == 0);
  REQUIRE(lattice.site(0).edge(1).offset == std::array<int, 2>{1, 0});
  REQUIRE_THAT(lattice.site(0).edge(1).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows() == 1);
  REQUIRE(adj_matrix.cols() == 1);
  REQUIRE(adj_matrix(0, 0) == 1);

  Matrix<std::complex<double>> hopping_matrix = lattice.HoppingMatrix(Vector<2>{0.5, 0.8});
  REQUIRE(hopping_matrix.rows() == 1);
  REQUIRE(hopping_matrix.cols() == 1);
  REQUIRE(hopping_matrix(0, 0) == std::complex<double>(2.0 * cos(0.5) + 2.0 * cos(0.8), 0));

  NumericArray<double> w(1);
  Matrix<std::complex<double>> v(1, 1);
  diagonalize_hermitian(hopping_matrix, w, v);
  REQUIRE_THAT(w[0], Catch::Matchers::WithinAbs(2.0 * cos(0.5) + 2.0 * cos(0.8), 1e-10));
}

class GrapheneLatticeTest : public TwoDimensionalLattice {
 public:
  GrapheneLatticeTest()
      : TwoDimensionalLattice(Vector<2>(1.5, 0.5 * sqrt(3.0)), Vector<2>(1.5, -0.5 * sqrt(3.0))) {
    add_site(Point<2>{0, 0});
    add_site(Point<2>{0.5, 0.5 * sqrt(3.0)});
    add_edge(0, 1, {0, 0}, 1.0);
    add_edge(1, 0, {1, 0}, 1.0);
    add_edge(1, 0, {1, -1}, 1.0);
  }
};

TEST_CASE("GrapheneLattice", "[lattice]") {
  GrapheneLatticeTest lattice;
  REQUIRE(lattice.a1() == Vector<2>(1.5, 0.5 * sqrt(3.0)));
  REQUIRE(lattice.a2() == Vector<2>(1.5, -0.5 * sqrt(3.0)));
  REQUIRE_THAT(lattice.b1(), ApproxEqualVec(Vector<2>(2.0 * M_PI / 3.0, 2.0 * M_PI / sqrt(3.0))));
  REQUIRE_THAT(lattice.b2(), ApproxEqualVec(Vector<2>(2.0 * M_PI / 3.0, -2.0 * M_PI / sqrt(3.0))));

  REQUIRE(lattice.size() == 2);

  REQUIRE(lattice.site(0).position == Point<2>{0, 0});
  REQUIRE(lattice.site(0).size() == 1);
  REQUIRE(lattice.site(0).edge(0).dst == 1);
  REQUIRE(lattice.site(0).edge(0).offset == std::array<int, 2>{0, 0});
  REQUIRE_THAT(lattice.site(0).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  REQUIRE(lattice.site(1).position == Point<2>{0.5, 0.5 * sqrt(3.0)});
  REQUIRE(lattice.site(1).size() == 2);
  REQUIRE(lattice.site(1).edge(0).dst == 0);
  REQUIRE(lattice.site(1).edge(0).offset == std::array<int, 2>{1, 0});
  REQUIRE_THAT(lattice.site(1).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(1).edge(1).dst == 0);
  REQUIRE(lattice.site(1).edge(1).offset == std::array<int, 2>{1, -1});

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows() == 2);
  REQUIRE(adj_matrix.cols() == 2);
  REQUIRE(adj_matrix(0, 0) == 0);
  REQUIRE(adj_matrix(0, 1) == 1);
  REQUIRE(adj_matrix(1, 0) == 1);
  REQUIRE(adj_matrix(1, 1) == 0);

  Vector<2> d1{0.5, 0.5 * sqrt(3.0)};
  Vector<2> d2{0.5, -0.5 * sqrt(3.0)};
  Vector<2> d3{-1.0, 0.0};
  Vector<2> k{0.5, 0.8};
  std::complex<double> comp = std::complex<double>(0.0, 1.0);

  Matrix<std::complex<double>> hopping_matrix = lattice.HoppingMatrix(k);
  REQUIRE(hopping_matrix.rows() == 2);
  REQUIRE(hopping_matrix.cols() == 2);
  REQUIRE(hopping_matrix(0, 0) == std::complex<double>(0, 0));
  REQUIRE(hopping_matrix(0, 1) ==
          std::exp(comp * k.dot(d1)) + std::exp(comp * k.dot(d2)) + std::exp(comp * k.dot(d3)));
  REQUIRE(hopping_matrix(1, 0) ==
          std::exp(-comp * k.dot(d1)) + std::exp(-comp * k.dot(d2)) + std::exp(-comp * k.dot(d3)));
  REQUIRE(hopping_matrix(1, 1) == std::complex<double>(0, 0));

  NumericArray<double> w(2);
  Matrix<std::complex<double>> v(2, 2);
  diagonalize_hermitian(hopping_matrix, w, v);
  double f_1 = 2.0 * cos(sqrt(3.0) * 0.8) + 4.0 * cos(0.5 * sqrt(3.0) * 0.8) * cos(1.5 * 0.5);
  REQUIRE_THAT(w[0], Catch::Matchers::WithinAbs(-sqrt(3.0 + f_1), 1e-10));
  REQUIRE_THAT(w[1], Catch::Matchers::WithinAbs(sqrt(3.0 + f_1), 1e-10));
}

class GrapheneLatticeExtendedTest : public TwoDimensionalLattice {
 public:
  GrapheneLatticeExtendedTest()
      : TwoDimensionalLattice(Vector<2>(3.0, 0), Vector<2>(0, sqrt(3.0))) {
    add_site(Point<2>{0, 0});
    add_site(Point<2>{0.5, 0.5 * sqrt(3.0)});
    add_site(Point<2>{1.5, 0.5 * sqrt(3.0)});
    add_site(Point<2>{2.0, 0});

    add_edge(0, 1, {0, 0}, 1.0);
    add_edge(1, 2, {0, 0}, 1.0);
    add_edge(2, 3, {0, 0}, 1.0);

    add_edge(3, 0, {1, 0}, 1.0);
    add_edge(1, 0, {0, 1}, 1.0);
    add_edge(2, 3, {0, 1}, 1.0);
  }
};

TEST_CASE("GrapheneLatticeExtended", "[lattice]") {
  GrapheneLatticeExtendedTest lattice;
  REQUIRE(lattice.size() == 4);

  REQUIRE(lattice.site(0).position == Point<2>{0, 0});
  REQUIRE(lattice.site(0).size() == 1);
  REQUIRE(lattice.site(0).edge(0).dst == 1);
  REQUIRE(lattice.site(0).edge(0).offset == std::array<int, 2>{0, 0});
  REQUIRE_THAT(lattice.site(0).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  REQUIRE(lattice.site(1).position == Point<2>{0.5, 0.5 * sqrt(3.0)});
  REQUIRE(lattice.site(1).size() == 2);
  REQUIRE(lattice.site(1).edge(0).dst == 2);
  REQUIRE(lattice.site(1).edge(0).offset == std::array<int, 2>{0, 0});
  REQUIRE_THAT(lattice.site(1).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(1).edge(1).dst == 0);
  REQUIRE(lattice.site(1).edge(1).offset == std::array<int, 2>{0, 1});
  REQUIRE_THAT(lattice.site(1).edge(1).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  REQUIRE(lattice.site(2).position == Point<2>{1.5, 0.5 * sqrt(3.0)});
  REQUIRE(lattice.site(2).size() == 2);
  REQUIRE(lattice.site(2).edge(0).dst == 3);
  REQUIRE(lattice.site(2).edge(0).offset == std::array<int, 2>{0, 0});
  REQUIRE_THAT(lattice.site(2).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(2).edge(1).dst == 3);
  REQUIRE(lattice.site(2).edge(1).offset == std::array<int, 2>{0, 1});
  REQUIRE_THAT(lattice.site(2).edge(1).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  REQUIRE(lattice.site(3).position == Point<2>{2.0, 0});
  REQUIRE(lattice.site(3).size() == 1);
  REQUIRE(lattice.site(3).edge(0).dst == 0);
  REQUIRE(lattice.site(3).edge(0).offset == std::array<int, 2>{1, 0});
  REQUIRE_THAT(lattice.site(3).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows() == 4);
  REQUIRE(adj_matrix.cols() == 4);
  REQUIRE(adj_matrix(0, 0) == 0);
  REQUIRE(adj_matrix(0, 1) == 1);
  REQUIRE(adj_matrix(0, 2) == 0);
  REQUIRE(adj_matrix(0, 3) == 1);
  REQUIRE(adj_matrix(1, 0) == 1);
  REQUIRE(adj_matrix(1, 1) == 0);
  REQUIRE(adj_matrix(1, 2) == 1);
  REQUIRE(adj_matrix(1, 3) == 0);
  REQUIRE(adj_matrix(2, 0) == 0);
  REQUIRE(adj_matrix(2, 1) == 1);
  REQUIRE(adj_matrix(2, 2) == 0);
  REQUIRE(adj_matrix(2, 3) == 1);
  REQUIRE(adj_matrix(3, 0) == 1);
  REQUIRE(adj_matrix(3, 1) == 0);
  REQUIRE(adj_matrix(3, 2) == 1);
  REQUIRE(adj_matrix(3, 3) == 0);
}

class TriangularLatticeTest : public TwoDimensionalLattice {
 public:
  TriangularLatticeTest(double a = 1.0)
      : TwoDimensionalLattice(Vector<2>(a, 0), Vector<2>(0.5 * a, 0.5 * a * sqrt(3.0))) {
    add_site(Point<2>{0, 0});
    add_edge(0, 0, {1, 0}, 1.0);
    add_edge(0, 0, {0, 1}, 1.0);
    add_edge(0, 0, {1, -1}, 1.0);
  }
};

TEST_CASE("TriangularLattice", "[lattice]") {
  TriangularLatticeTest lattice;
  REQUIRE(lattice.size() == 1);

  REQUIRE(lattice.site(0).position == Point<2>{0, 0});
  REQUIRE(lattice.site(0).size() == 3);
  REQUIRE(lattice.site(0).edge(0).dst == 0);
  REQUIRE(lattice.site(0).edge(0).offset == std::array<int, 2>{1, 0});
  REQUIRE_THAT(lattice.site(0).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(0).edge(1).dst == 0);
  REQUIRE(lattice.site(0).edge(1).offset == std::array<int, 2>{0, 1});
  REQUIRE_THAT(lattice.site(0).edge(1).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(0).edge(2).dst == 0);
  REQUIRE(lattice.site(0).edge(2).offset == std::array<int, 2>{1, -1});
  REQUIRE_THAT(lattice.site(0).edge(2).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows() == 1);
  REQUIRE(adj_matrix.cols() == 1);
  REQUIRE(adj_matrix(0, 0) == 1);
}

class KagomeLatticeTest : public TwoDimensionalLattice {
 public:
  KagomeLatticeTest() : TwoDimensionalLattice(Vector<2>(2, 0), Vector<2>(1.0, sqrt(3.0))) {
    add_site(Point<2>{0, 0});
    add_site(Point<2>{1.0, 0});
    add_site(Point<2>{0.5, 0.5 * sqrt(3.0)});
    add_edge(0, 1, {0, 0}, 1.0);
    add_edge(1, 2, {0, 0}, 1.0);
    add_edge(2, 0, {0, 0}, 1.0);
    add_edge(1, 0, {1, 0}, 1.0);
    add_edge(2, 0, {0, 1}, 1.0);
    add_edge(1, 2, {1, -1}, 1.0);
  }
};

TEST_CASE("KagomeLattice", "[lattice]") {
  KagomeLatticeTest lattice;
  REQUIRE(lattice.size() == 3);

  REQUIRE(lattice.site(0).position == Point<2>{0, 0});
  REQUIRE(lattice.site(0).size() == 1);
  REQUIRE(lattice.site(0).edge(0).dst == 1);
  REQUIRE(lattice.site(0).edge(0).offset == std::array<int, 2>{0, 0});
  REQUIRE_THAT(lattice.site(0).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  REQUIRE(lattice.site(1).position == Point<2>{1.0, 0});
  REQUIRE(lattice.site(1).size() == 3);
  REQUIRE(lattice.site(1).edge(0).dst == 2);
  REQUIRE(lattice.site(1).edge(0).offset == std::array<int, 2>{0, 0});
  REQUIRE_THAT(lattice.site(1).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(1).edge(1).dst == 0);
  REQUIRE(lattice.site(1).edge(1).offset == std::array<int, 2>{1, 0});
  REQUIRE_THAT(lattice.site(1).edge(1).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(1).edge(2).dst == 2);
  REQUIRE(lattice.site(1).edge(2).offset == std::array<int, 2>{1, -1});
  REQUIRE_THAT(lattice.site(1).edge(2).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  REQUIRE(lattice.site(2).position == Point<2>{0.5, 0.5 * sqrt(3.0)});
  REQUIRE(lattice.site(2).size() == 2);
  REQUIRE(lattice.site(2).edge(0).dst == 0);
  REQUIRE(lattice.site(2).edge(0).offset == std::array<int, 2>{0, 0});
  REQUIRE_THAT(lattice.site(2).edge(0).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE(lattice.site(2).edge(1).dst == 0);
  REQUIRE(lattice.site(2).edge(1).offset == std::array<int, 2>{0, 1});
  REQUIRE_THAT(lattice.site(2).edge(1).weight, Catch::Matchers::WithinAbs(1.0, 1e-10));

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows() == 3);
  REQUIRE(adj_matrix.cols() == 3);
  REQUIRE(adj_matrix(0, 0) == 0);
  REQUIRE(adj_matrix(0, 1) == 1);
  REQUIRE(adj_matrix(0, 2) == 1);
  REQUIRE(adj_matrix(1, 0) == 1);
  REQUIRE(adj_matrix(1, 1) == 0);
  REQUIRE(adj_matrix(1, 2) == 1);
  REQUIRE(adj_matrix(2, 0) == 1);
  REQUIRE(adj_matrix(2, 1) == 1);
  REQUIRE(adj_matrix(2, 2) == 0);
}
