#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <complex>

#include "cctb/lattice.h"
#include "cctb/vec.h"

struct ApproxEqualVec : Catch::Matchers::MatcherGenericBase {
 public:
  ApproxEqualVec(Vec<float> v) : v_(v) {}

  bool match(Vec<float> const &v) const {
    if (v.size() != v_.size()) {
      return false;
    }
    for (std::size_t i = 0; i < v.size(); i++) {
      if (std::abs(v[i] - v_[i]) > 1e-6) {
        return false;
      }
    }
    return true;
  }

  std::string describe() const override { return "ApproxEqualVec"; }

 private:
  Vec<float> v_;
};

class LinearChainTest : public OneDimensionalLattice {
 public:
  LinearChainTest(int size) : OneDimensionalLattice(Vec<float>{1.0f * size}) {
    for (int i = 0; i < size; i++) {
      AddSite(Site(Vec<float>{static_cast<float>(i)}));
    }
    for (int i = 0; i < size - 1; i++) {
      AddEdge(Edge({0}, i, i + 1));
    }
    // if periodic, add the last edge
    AddEdge(Edge({1}, size - 1, 0));
  }
};

TEST_CASE("LinearChain", "[lattice]") {
  LinearChainTest lattice(2);
  REQUIRE(lattice.Size() == 2);
  REQUIRE(lattice.SiteAt(0).position == Vec<float>{0});
  REQUIRE(lattice.SiteAt(1).position == Vec<float>{1});
  REQUIRE(lattice.Edges().size() == 2);
  REQUIRE(lattice.Edges()[0].relative_index == std::vector<int>{0});
  REQUIRE(lattice.Edges()[1].relative_index == std::vector<int>{1});
  REQUIRE(lattice.Edges()[0].from_index == 0);
  REQUIRE(lattice.Edges()[0].to_index == 1);
  REQUIRE(lattice.Edges()[1].from_index == 1);
  REQUIRE(lattice.Edges()[1].to_index == 0);

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows == 2);
  REQUIRE(adj_matrix.cols == 2);
  REQUIRE(adj_matrix(0, 0) == 0);
  REQUIRE(adj_matrix(0, 1) == 1);
  REQUIRE(adj_matrix(1, 0) == 1);
  REQUIRE(adj_matrix(1, 1) == 0);

  Matrix<std::complex<float>> hopping_matrix = lattice.HoppingMatrix(0.5);
  REQUIRE(hopping_matrix.rows == 2);
  REQUIRE(hopping_matrix.cols == 2);
  REQUIRE(hopping_matrix(0, 0) == std::complex<float>(0, 0));
  REQUIRE(hopping_matrix(0, 1) == std::complex<float>(2.0f * cosf(0.5), 0));
  REQUIRE(hopping_matrix(1, 0) == std::complex<float>(2.0f * cosf(0.5), 0));
  REQUIRE(hopping_matrix(1, 1) == std::complex<float>(0, 0));
}

class SquareLatticeTest : public TwoDimensionalLattice {
 public:
  SquareLatticeTest() : TwoDimensionalLattice(Vec<float>{1.0f, 0}, Vec<float>{0, 1.0f}) {
    AddSite(Site(Vec<float>{0, 0}));
    AddEdge(Edge({0, 1}, 0, 0));
    AddEdge(Edge({1, 0}, 0, 0));
  }
};

TEST_CASE("SquareLattice", "[lattice]") {
  SquareLatticeTest lattice;
  REQUIRE(lattice.Size() == 1);
  REQUIRE(lattice.SiteAt(0).position == Vec<float>{0, 0});
  REQUIRE(lattice.Edges().size() == 2);
  REQUIRE(lattice.Edges()[0].relative_index == std::vector<int>{0, 1});
  REQUIRE(lattice.Edges()[1].relative_index == std::vector<int>{1, 0});
  REQUIRE(lattice.Edges()[0].from_index == 0);
  REQUIRE(lattice.Edges()[0].to_index == 0);
  REQUIRE(lattice.Edges()[1].from_index == 0);
  REQUIRE(lattice.Edges()[1].to_index == 0);

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows == 1);
  REQUIRE(adj_matrix.cols == 1);
  REQUIRE(adj_matrix(0, 0) == 1);

  Matrix<std::complex<float>> hopping_matrix = lattice.HoppingMatrix(Vec<float>{0.5, 0.8});
  REQUIRE(hopping_matrix.rows == 1);
  REQUIRE(hopping_matrix.cols == 1);
  REQUIRE(hopping_matrix(0, 0) == std::complex<float>(2.0f * cosf(0.5) + 2.0f * cosf(0.8), 0));
}

class GrapheneLatticeTest : public TwoDimensionalLattice {
 public:
  GrapheneLatticeTest()
      : TwoDimensionalLattice(Vec<float>(1.5f, 0.5f * sqrtf(3.0f)),
                              Vec<float>(1.5f, -0.5f * sqrtf(3.0f))) {
    AddSite(Site(Vec<float>{0, 0}));
    AddSite(Site(Vec<float>{0.5, 0.5f * sqrtf(3.0f)}));
    AddEdge(Edge({0, 0}, 0, 1));
    AddEdge(Edge({1, 0}, 1, 0));
    AddEdge(Edge({1, -1}, 1, 0));
  }
};

TEST_CASE("GrapheneLattice", "[lattice]") {
  GrapheneLatticeTest lattice;
  REQUIRE(lattice.Size() == 2);
  REQUIRE(lattice.a1() == Vec<float>(1.5f, 0.5f * sqrtf(3.0f)));
  REQUIRE(lattice.a2() == Vec<float>(1.5f, -0.5f * sqrtf(3.0f)));
  REQUIRE_THAT(lattice.b1(),
               ApproxEqualVec(Vec<float>(2.0f * M_PI / 3.0f, 2.0f * M_PI / sqrtf(3.0f))));
  REQUIRE_THAT(lattice.b2(),
               ApproxEqualVec(Vec<float>(2.0f * M_PI / 3.0f, -2.0f * M_PI / sqrtf(3.0f))));

  REQUIRE(lattice.SiteAt(0).position == Vec<float>{0, 0});
  REQUIRE(lattice.SiteAt(1).position == Vec<float>{0.5, 0.5f * sqrtf(3.0f)});
  REQUIRE(lattice.Edges().size() == 3);
  REQUIRE(lattice.Edges()[0].relative_index == std::vector<int>{0, 0});
  REQUIRE(lattice.Edges()[1].relative_index == std::vector<int>{1, 0});
  REQUIRE(lattice.Edges()[2].relative_index == std::vector<int>{1, -1});
  REQUIRE(lattice.Edges()[0].from_index == 0);
  REQUIRE(lattice.Edges()[0].to_index == 1);
  REQUIRE(lattice.Edges()[1].from_index == 1);
  REQUIRE(lattice.Edges()[1].to_index == 0);
  REQUIRE(lattice.Edges()[2].from_index == 1);
  REQUIRE(lattice.Edges()[2].to_index == 0);

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows == 2);
  REQUIRE(adj_matrix.cols == 2);
  REQUIRE(adj_matrix(0, 0) == 0);
  REQUIRE(adj_matrix(0, 1) == 1);
  REQUIRE(adj_matrix(1, 0) == 1);
  REQUIRE(adj_matrix(1, 1) == 0);

  Vec<float> d1{0.5f, 0.5f * sqrtf(3.0)};
  Vec<float> d2{0.5f, -0.5f * sqrtf(3.0)};
  Vec<float> d3{-1.0f, 0.0f};
  Vec<float> k{0.5, 0.8};
  std::complex<float> comp = std::complex<float>(0.0f, 1.0f);

  Matrix<std::complex<float>> hopping_matrix = lattice.HoppingMatrix(k);
  REQUIRE(hopping_matrix.rows == 2);
  REQUIRE(hopping_matrix.cols == 2);
  REQUIRE(hopping_matrix(0, 0) == std::complex<float>(0, 0));
  REQUIRE(hopping_matrix(0, 1) ==
          std::exp(comp * k.dot(d1)) + std::exp(comp * k.dot(d2)) + std::exp(comp * k.dot(d3)));
  REQUIRE(hopping_matrix(1, 0) ==
          std::exp(-comp * k.dot(d1)) + std::exp(-comp * k.dot(d2)) + std::exp(-comp * k.dot(d3)));
  REQUIRE(hopping_matrix(1, 1) == std::complex<float>(0, 0));
}

class GrapheneLatticeExtendedTest : public TwoDimensionalLattice {
 public:
  GrapheneLatticeExtendedTest()
      : TwoDimensionalLattice(Vec<float>(3.0f, 0), Vec<float>(0, sqrtf(3.0f))) {
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

TEST_CASE("GrapheneLatticeExtended", "[lattice]") {
  GrapheneLatticeExtendedTest lattice;
  REQUIRE(lattice.Size() == 4);
  REQUIRE(lattice.SiteAt(0).position == Vec<float>{0, 0});
  REQUIRE(lattice.SiteAt(1).position == Vec<float>{0.5f, 0.5f * sqrtf(3.0f)});
  REQUIRE(lattice.SiteAt(2).position == Vec<float>{1.5f, 0.5f * sqrtf(3.0f)});
  REQUIRE(lattice.SiteAt(3).position == Vec<float>{2.0f, 0});
  REQUIRE(lattice.Edges().size() == 6);
  REQUIRE(lattice.Edges()[0].relative_index == std::vector<int>{0, 0});
  REQUIRE(lattice.Edges()[1].relative_index == std::vector<int>{0, 0});
  REQUIRE(lattice.Edges()[2].relative_index == std::vector<int>{0, 0});
  REQUIRE(lattice.Edges()[3].relative_index == std::vector<int>{1, 0});
  REQUIRE(lattice.Edges()[4].relative_index == std::vector<int>{0, 1});
  REQUIRE(lattice.Edges()[5].relative_index == std::vector<int>{0, 1});
  REQUIRE(lattice.Edges()[0].from_index == 0);
  REQUIRE(lattice.Edges()[0].to_index == 1);
  REQUIRE(lattice.Edges()[1].from_index == 1);
  REQUIRE(lattice.Edges()[1].to_index == 2);
  REQUIRE(lattice.Edges()[2].from_index == 2);
  REQUIRE(lattice.Edges()[2].to_index == 3);
  REQUIRE(lattice.Edges()[3].from_index == 3);
  REQUIRE(lattice.Edges()[3].to_index == 0);
  REQUIRE(lattice.Edges()[4].from_index == 1);
  REQUIRE(lattice.Edges()[4].to_index == 0);
  REQUIRE(lattice.Edges()[5].from_index == 2);
  REQUIRE(lattice.Edges()[5].to_index == 3);

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows == 4);
  REQUIRE(adj_matrix.cols == 4);
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
  TriangularLatticeTest(float a = 1.0f)
      : TwoDimensionalLattice(Vec<float>(a, 0), Vec<float>(0.5f * a, 0.5f * a * sqrtf(3.0f))) {
    AddSite(Site(Vec<float>{0, 0}));
    AddEdge(Edge({1, 0}, 0, 0));
    AddEdge(Edge({0, 1}, 0, 0));
    AddEdge(Edge({1, -1}, 0, 0));
  }
};

TEST_CASE("TriangularLattice", "[lattice]") {
  TriangularLatticeTest lattice;
  REQUIRE(lattice.Size() == 1);
  REQUIRE(lattice.SiteAt(0).position == Vec<float>{0, 0});
  REQUIRE(lattice.Edges().size() == 3);
  REQUIRE(lattice.Edges()[0].relative_index == std::vector<int>{1, 0});
  REQUIRE(lattice.Edges()[1].relative_index == std::vector<int>{0, 1});
  REQUIRE(lattice.Edges()[2].relative_index == std::vector<int>{1, -1});
  REQUIRE(lattice.Edges()[0].from_index == 0);
  REQUIRE(lattice.Edges()[0].to_index == 0);
  REQUIRE(lattice.Edges()[1].from_index == 0);
  REQUIRE(lattice.Edges()[1].to_index == 0);
  REQUIRE(lattice.Edges()[2].from_index == 0);
  REQUIRE(lattice.Edges()[2].to_index == 0);

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows == 1);
  REQUIRE(adj_matrix.cols == 1);
  REQUIRE(adj_matrix(0, 0) == 1);
}

class KagomeLatticeTest : public TwoDimensionalLattice {
 public:
  KagomeLatticeTest() : TwoDimensionalLattice(Vec<float>(2, 0), Vec<float>(1.0f, sqrtf(3.0f))) {
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

TEST_CASE("KagomeLattice", "[lattice]") {
  KagomeLatticeTest lattice;
  REQUIRE(lattice.Size() == 3);
  REQUIRE(lattice.SiteAt(0).position == Vec<float>{0, 0});
  REQUIRE(lattice.SiteAt(1).position == Vec<float>{1.0f, 0});
  REQUIRE(lattice.SiteAt(2).position == Vec<float>{0.5f, 0.5f * sqrtf(3.0f)});
  REQUIRE(lattice.Edges().size() == 6);
  REQUIRE(lattice.Edges()[0].relative_index == std::vector<int>{0, 0});
  REQUIRE(lattice.Edges()[1].relative_index == std::vector<int>{0, 0});
  REQUIRE(lattice.Edges()[2].relative_index == std::vector<int>{0, 0});
  REQUIRE(lattice.Edges()[3].relative_index == std::vector<int>{1, 0});
  REQUIRE(lattice.Edges()[4].relative_index == std::vector<int>{0, 1});
  REQUIRE(lattice.Edges()[5].relative_index == std::vector<int>{1, -1});
  REQUIRE(lattice.Edges()[0].from_index == 0);
  REQUIRE(lattice.Edges()[0].to_index == 1);
  REQUIRE(lattice.Edges()[1].from_index == 1);
  REQUIRE(lattice.Edges()[1].to_index == 2);
  REQUIRE(lattice.Edges()[2].from_index == 2);
  REQUIRE(lattice.Edges()[2].to_index == 0);
  REQUIRE(lattice.Edges()[3].from_index == 1);
  REQUIRE(lattice.Edges()[3].to_index == 0);
  REQUIRE(lattice.Edges()[4].from_index == 2);
  REQUIRE(lattice.Edges()[4].to_index == 0);
  REQUIRE(lattice.Edges()[5].from_index == 1);
  REQUIRE(lattice.Edges()[5].to_index == 2);

  Matrix<int> adj_matrix = lattice.AdjMatrix();
  REQUIRE(adj_matrix.rows == 3);
  REQUIRE(adj_matrix.cols == 3);
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