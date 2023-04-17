#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <complex>

#include "Lattice/Lattice.h"
#include "TestMatchers.h"

class LinearChainTest : public OneDimensionalLattice<1> {
 public:
  using NumberType = double;
  using ComplexType = std::complex<NumberType>;
  using VectorType = Eigen::Vector<NumberType, 1>;
  using MatrixType = Eigen::Matrix<NumberType, 1, 1>;

  LinearChainTest(size_t size)
      : OneDimensionalLattice(VectorType(1.0 * static_cast<double>(size))) {
    MatrixType onsite = MatrixType::Zero();
    MatrixType hopping = MatrixType::Identity();

    for (size_t i = 0; i < size; i++) {
      add_site(i, VectorType(static_cast<double>(i)), onsite);
    }
    for (size_t i = 0; i < size - 1; i++) {
      add_edge(i, i + 1, {0}, hopping);
    }
    // if periodic, add the last edge
    add_edge(size - 1, 0, {1}, hopping);
  }
};

TEST_CASE("LinearChain", "[lattice]") {
  LinearChainTest lattice(2);
  REQUIRE(lattice.size() == 2);
  REQUIRE(lattice.site(0).index() == 0);
  REQUIRE(lattice.site(0).position() == LinearChainTest::VectorType(0));
  REQUIRE(lattice.edge(0).src() == 0);
  REQUIRE(lattice.edge(0).dst() == 1);
  REQUIRE(lattice.edge(0).basis_index(0) == 0);
  REQUIRE(lattice.site(1).index() == 1);
  REQUIRE(lattice.site(1).position() == LinearChainTest::VectorType(1));
  REQUIRE(lattice.edge(1).src() == 1);
  REQUIRE(lattice.edge(1).dst() == 0);
  REQUIRE(lattice.edge(1).basis_index(0) == 1);

  Eigen::MatrixXi adj_matrix = lattice.graph().adj_matrix();
  REQUIRE(adj_matrix.rows() == 2);
  REQUIRE(adj_matrix.cols() == 2);
  REQUIRE(adj_matrix(0, 0) == 0);
  REQUIRE(adj_matrix(0, 1) == 1);
  REQUIRE(adj_matrix(1, 0) == 1);
  REQUIRE(adj_matrix(1, 1) == 0);

  Eigen::MatrixXcd hopping_matrix = lattice.hopping_matrix(LinearChainTest::VectorType(0.5));
  REQUIRE(hopping_matrix.rows() == 2);
  REQUIRE(hopping_matrix.cols() == 2);
  REQUIRE(hopping_matrix(0, 0) == std::complex<double>(0, 0));
  REQUIRE(hopping_matrix(0, 1) == std::complex<double>(2.0 * cos(0.5), 0));
  REQUIRE(hopping_matrix(1, 0) == std::complex<double>(2.0 * cos(0.5), 0));
  REQUIRE(hopping_matrix(1, 1) == std::complex<double>(0, 0));

  Eigen::VectorXd w(2);
  Eigen::MatrixXcd v(2, 2);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(hopping_matrix);
  w = es.eigenvalues();
  v = es.eigenvectors();

  REQUIRE_THAT(w[0], Catch::Matchers::WithinAbs(-2.0 * cos(0.5), 1e-10));
  REQUIRE_THAT(w[1], Catch::Matchers::WithinAbs(2.0 * cos(0.5), 1e-10));
}

class SquareLatticeTest : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using ComplexType = std::complex<NumberType>;
  using VectorType = Eigen::Vector<NumberType, 2>;
  using MatrixType = Eigen::Matrix<NumberType, 1, 1>;

  SquareLatticeTest() : TwoDimensionalLattice(VectorType(1.0, 0), VectorType(0, 1.0)) {
    MatrixType onsite = MatrixType::Zero();
    MatrixType hopping = MatrixType::Identity();

    add_site(0, VectorType(0, 0), onsite);
    add_edge(0, 0, {0, 1}, hopping);
    add_edge(0, 0, {1, 0}, hopping);
  }
};

TEST_CASE("SquareLattice", "[lattice]") {
  SquareLatticeTest lattice;
  REQUIRE(lattice.size() == 1);

  REQUIRE(lattice.site(0).index() == 0);
  REQUIRE(lattice.site(0).position() == SquareLatticeTest::VectorType(0, 0));
  REQUIRE(lattice.edge(0).src() == 0);
  REQUIRE(lattice.edge(0).dst() == 0);
  REQUIRE(lattice.edge(0).basis_index(0) == 0);
  REQUIRE(lattice.edge(0).basis_index(1) == 1);
  REQUIRE(lattice.edge(1).src() == 0);
  REQUIRE(lattice.edge(1).dst() == 0);
  REQUIRE(lattice.edge(1).basis_index(0) == 1);
  REQUIRE(lattice.edge(1).basis_index(1) == 0);

  Eigen::MatrixXi adj_matrix = lattice.graph().adj_matrix();
  REQUIRE(adj_matrix.rows() == 1);
  REQUIRE(adj_matrix.cols() == 1);
  REQUIRE(adj_matrix(0, 0) == 1);

  Eigen::MatrixXcd hopping_matrix = lattice.hopping_matrix(SquareLatticeTest::VectorType(0.5, 0.8));
  REQUIRE(hopping_matrix.rows() == 1);
  REQUIRE(hopping_matrix.cols() == 1);
  REQUIRE(hopping_matrix(0, 0) == std::complex<double>(2.0 * cos(0.5) + 2.0 * cos(0.8), 0));

  Eigen::VectorXd w(1);
  Eigen::MatrixXcd v(1, 1);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(hopping_matrix);
  w = es.eigenvalues();
  v = es.eigenvectors();
  REQUIRE_THAT(w[0], Catch::Matchers::WithinAbs(2.0 * cos(0.5) + 2.0 * cos(0.8), 1e-10));
}

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

TEST_CASE("GrapheneLattice", "[lattice]") {
  GrapheneLatticeTest lattice;
  REQUIRE(lattice.lattice_vector(0) == GrapheneLatticeTest::VectorType(1.5, 0.5 * sqrt(3.0)));
  REQUIRE(lattice.lattice_vector(1) == GrapheneLatticeTest::VectorType(1.5, -0.5 * sqrt(3.0)));
  REQUIRE_THAT((lattice.reciprocal_vector(0) -
                GrapheneLatticeTest::VectorType(2.0 * M_PI / 3.0, 2.0 * M_PI / sqrt(3.0)))
                   .norm(),
               Catch::Matchers::WithinAbs(0.0, 1e-10));
  REQUIRE_THAT((lattice.reciprocal_vector(1) -
                GrapheneLatticeTest::VectorType(2.0 * M_PI / 3.0, -2.0 * M_PI / sqrt(3.0)))
                   .norm(),
               Catch::Matchers::WithinAbs(0.0, 1e-10));

  REQUIRE(lattice.size() == 2);
  REQUIRE(lattice.site(0).index() == 0);
  REQUIRE(lattice.site(0).position() == GrapheneLatticeTest::VectorType(0, 0));
  REQUIRE(lattice.site(1).index() == 1);
  REQUIRE(lattice.site(1).position() == GrapheneLatticeTest::VectorType(0.5, 0.5 * sqrt(3.0)));

  REQUIRE(lattice.edge(0).src() == 0);
  REQUIRE(lattice.edge(0).dst() == 1);
  REQUIRE(lattice.edge(0).basis_index(0) == 0);
  REQUIRE(lattice.edge(0).basis_index(1) == 0);

  REQUIRE(lattice.edge(1).src() == 1);
  REQUIRE(lattice.edge(1).dst() == 0);
  REQUIRE(lattice.edge(1).basis_index(0) == 1);
  REQUIRE(lattice.edge(1).basis_index(1) == 0);

  REQUIRE(lattice.edge(2).src() == 1);
  REQUIRE(lattice.edge(2).dst() == 0);
  REQUIRE(lattice.edge(2).basis_index(0) == 1);
  REQUIRE(lattice.edge(2).basis_index(1) == -1);

  Eigen::MatrixXi adj_matrix = lattice.graph().adj_matrix();
  REQUIRE(adj_matrix.rows() == 2);
  REQUIRE(adj_matrix.cols() == 2);
  REQUIRE(adj_matrix(0, 0) == 0);
  REQUIRE(adj_matrix(0, 1) == 1);
  REQUIRE(adj_matrix(1, 0) == 1);
  REQUIRE(adj_matrix(1, 1) == 0);

  GrapheneLatticeTest::VectorType d1(0.5, 0.5 * sqrt(3.0));
  GrapheneLatticeTest::VectorType d2(0.5, -0.5 * sqrt(3.0));
  GrapheneLatticeTest::VectorType d3(-1.0, 0.0);
  GrapheneLatticeTest::VectorType k(0.5, 0.8);
  std::complex<double> comp = std::complex<double>(0.0, 1.0);

  Eigen::MatrixXcd hopping_matrix = lattice.hopping_matrix(k);
  REQUIRE(hopping_matrix.rows() == 2);
  REQUIRE(hopping_matrix.cols() == 2);
  REQUIRE(hopping_matrix(0, 0) == std::complex<double>(0, 0));
  REQUIRE(hopping_matrix(0, 1) ==
          std::exp(comp * k.dot(d1)) + std::exp(comp * k.dot(d2)) + std::exp(comp * k.dot(d3)));
  REQUIRE(hopping_matrix(1, 0) ==
          std::exp(-comp * k.dot(d1)) + std::exp(-comp * k.dot(d2)) + std::exp(-comp * k.dot(d3)));
  REQUIRE(hopping_matrix(1, 1) == std::complex<double>(0, 0));

  Eigen::VectorXd w(2);
  Eigen::MatrixXcd v(2, 2);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(hopping_matrix);
  w = es.eigenvalues();
  v = es.eigenvectors();
  double f_1 = 2.0 * cos(sqrt(3.0) * 0.8) + 4.0 * cos(0.5 * sqrt(3.0) * 0.8) * cos(1.5 * 0.5);
  REQUIRE_THAT(w[0], Catch::Matchers::WithinAbs(-sqrt(3.0 + f_1), 1e-10));
  REQUIRE_THAT(w[1], Catch::Matchers::WithinAbs(sqrt(3.0 + f_1), 1e-10));
}

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

TEST_CASE("GrapheneLatticeExtended", "[lattice]") {
  GrapheneLatticeExtendedTest lattice;
  REQUIRE(lattice.size() == 4);

  REQUIRE(lattice.site(0).index() == 0);
  REQUIRE(lattice.site(0).position() == GrapheneLatticeExtendedTest::VectorType{0, 0});

  REQUIRE(lattice.site(1).index() == 1);
  REQUIRE(lattice.site(1).position() ==
          GrapheneLatticeExtendedTest::VectorType{0.5, 0.5 * sqrt(3.0)});

  REQUIRE(lattice.site(2).index() == 2);
  REQUIRE(lattice.site(2).position() ==
          GrapheneLatticeExtendedTest::VectorType{1.5, 0.5 * sqrt(3.0)});

  REQUIRE(lattice.site(3).index() == 3);
  REQUIRE(lattice.site(3).position() == GrapheneLatticeExtendedTest::VectorType{2.0, 0});

  REQUIRE(lattice.edge(0).src() == 0);
  REQUIRE(lattice.edge(0).dst() == 1);
  REQUIRE(lattice.edge(0).basis_index(0) == 0);
  REQUIRE(lattice.edge(0).basis_index(1) == 0);

  REQUIRE(lattice.edge(1).src() == 1);
  REQUIRE(lattice.edge(1).dst() == 2);
  REQUIRE(lattice.edge(1).basis_index(0) == 0);
  REQUIRE(lattice.edge(1).basis_index(1) == 0);

  REQUIRE(lattice.edge(2).src() == 2);
  REQUIRE(lattice.edge(2).dst() == 3);
  REQUIRE(lattice.edge(2).basis_index(0) == 0);
  REQUIRE(lattice.edge(2).basis_index(1) == 0);

  REQUIRE(lattice.edge(3).src() == 3);
  REQUIRE(lattice.edge(3).dst() == 0);
  REQUIRE(lattice.edge(3).basis_index(0) == 1);
  REQUIRE(lattice.edge(3).basis_index(1) == 0);

  REQUIRE(lattice.edge(4).src() == 1);
  REQUIRE(lattice.edge(4).dst() == 0);
  REQUIRE(lattice.edge(4).basis_index(0) == 0);
  REQUIRE(lattice.edge(4).basis_index(1) == 1);

  REQUIRE(lattice.edge(5).src() == 2);
  REQUIRE(lattice.edge(5).dst() == 3);
  REQUIRE(lattice.edge(5).basis_index(0) == 0);
  REQUIRE(lattice.edge(5).basis_index(1) == 1);

  Eigen::MatrixXi adj_matrix = lattice.graph().adj_matrix();
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

TEST_CASE("TriangularLattice", "[lattice]") {
  TriangularLatticeTest lattice;
  REQUIRE(lattice.size() == 1);

  REQUIRE(lattice.site(0).index() == 0);
  REQUIRE(lattice.site(0).position() == TriangularLatticeTest::VectorType{0, 0});

  REQUIRE(lattice.edge(0).src() == 0);
  REQUIRE(lattice.edge(0).dst() == 0);
  REQUIRE(lattice.edge(0).basis_index(0) == 1);
  REQUIRE(lattice.edge(0).basis_index(1) == 0);

  REQUIRE(lattice.edge(1).src() == 0);
  REQUIRE(lattice.edge(1).dst() == 0);
  REQUIRE(lattice.edge(1).basis_index(0) == 0);
  REQUIRE(lattice.edge(1).basis_index(1) == 1);

  REQUIRE(lattice.edge(2).src() == 0);
  REQUIRE(lattice.edge(2).dst() == 0);
  REQUIRE(lattice.edge(2).basis_index(0) == 1);
  REQUIRE(lattice.edge(2).basis_index(1) == -1);

  Eigen::MatrixXi adj_matrix = lattice.graph().adj_matrix();
  REQUIRE(adj_matrix.rows() == 1);
  REQUIRE(adj_matrix.cols() == 1);
  REQUIRE(adj_matrix(0, 0) == 1);

  // TODO: test eigenvalues with closed form
}

class KagomeLatticeTest : public TwoDimensionalLattice<1> {
 public:
  using NumberType = double;
  using ComplexType = std::complex<NumberType>;
  using VectorType = Eigen::Vector<NumberType, 2>;
  using MatrixType = Eigen::Matrix<NumberType, 1, 1>;

  KagomeLatticeTest() : TwoDimensionalLattice(VectorType(2, 0), VectorType(1.0, sqrt(3.0))) {
    MatrixType onsite = MatrixType::Zero();
    MatrixType hopping = MatrixType::Identity();
    add_site(0, VectorType(0, 0), onsite);
    add_site(1, VectorType(1.0, 0), onsite);
    add_site(2, VectorType(0.5, 0.5 * sqrt(3.0)), onsite);
    add_edge(0, 1, {0, 0}, hopping);
    add_edge(1, 2, {0, 0}, hopping);
    add_edge(2, 0, {0, 0}, hopping);
    add_edge(1, 0, {1, 0}, hopping);
    add_edge(2, 0, {0, 1}, hopping);
    add_edge(1, 2, {1, -1}, hopping);
  }
};

TEST_CASE("KagomeLattice", "[lattice]") {
  KagomeLatticeTest lattice;
  REQUIRE(lattice.size() == 3);

  REQUIRE(lattice.site(0).index() == 0);
  REQUIRE(lattice.site(0).position() == KagomeLatticeTest::VectorType{0, 0});

  REQUIRE(lattice.site(1).index() == 1);
  REQUIRE(lattice.site(1).position() == KagomeLatticeTest::VectorType{1, 0});

  REQUIRE(lattice.site(2).index() == 2);
  REQUIRE(lattice.site(2).position() == KagomeLatticeTest::VectorType{0.5, 0.5 * sqrt(3.0)});

  REQUIRE(lattice.edge(0).src() == 0);
  REQUIRE(lattice.edge(0).dst() == 1);
  REQUIRE(lattice.edge(0).basis_index(0) == 0);
  REQUIRE(lattice.edge(0).basis_index(1) == 0);

  REQUIRE(lattice.edge(1).src() == 1);
  REQUIRE(lattice.edge(1).dst() == 2);
  REQUIRE(lattice.edge(1).basis_index(0) == 0);
  REQUIRE(lattice.edge(1).basis_index(1) == 0);

  REQUIRE(lattice.edge(2).src() == 2);
  REQUIRE(lattice.edge(2).dst() == 0);
  REQUIRE(lattice.edge(2).basis_index(0) == 0);
  REQUIRE(lattice.edge(2).basis_index(1) == 0);

  REQUIRE(lattice.edge(3).src() == 1);
  REQUIRE(lattice.edge(3).dst() == 0);
  REQUIRE(lattice.edge(3).basis_index(0) == 1);
  REQUIRE(lattice.edge(3).basis_index(1) == 0);

  REQUIRE(lattice.edge(4).src() == 2);
  REQUIRE(lattice.edge(4).dst() == 0);
  REQUIRE(lattice.edge(4).basis_index(0) == 0);
  REQUIRE(lattice.edge(4).basis_index(1) == 1);

  REQUIRE(lattice.edge(5).src() == 1);
  REQUIRE(lattice.edge(5).dst() == 2);
  REQUIRE(lattice.edge(5).basis_index(0) == 1);
  REQUIRE(lattice.edge(5).basis_index(1) == -1);

  Eigen::MatrixXi adj_matrix = lattice.graph().adj_matrix();
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

  // TODO: test eigenvalues with closed form
}
