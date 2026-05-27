#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/Grid.hpp"

using namespace ffroom;
using Catch::Matchers::WithinAbs;

TEST_CASE("Grid construction", "[grid]") {
    Grid g(4, 5, 6, 4.0, 5.0, 6.0);
    REQUIRE(g.Nx == 4);
    REQUIRE(g.Ny == 5);
    REQUIRE(g.Nz == 6);
    REQUIRE_THAT(g.dx, WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(g.dy, WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(g.dz, WithinAbs(1.0, 1e-12));

    REQUIRE(g.u.size() == size_t(5*5*6));
    REQUIRE(g.v.size() == size_t(4*6*6));
    REQUIRE(g.w.size() == size_t(4*5*7));
    REQUIRE(g.p.size() == size_t(4*5*6));
}

TEST_CASE("Grid index round-trip", "[grid]") {
    Grid g(3, 4, 5, 1.0, 1.0, 1.0);
    // c_idx must be unique for all (i,j,k)
    std::set<int> seen;
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
    for (int k = 0; k < 5; k++) {
        int idx = g.c_idx(i,j,k);
        REQUIRE(idx >= 0);
        REQUIRE(idx < 3*4*5);
        REQUIRE(seen.find(idx) == seen.end());
        seen.insert(idx);
    }
}

TEST_CASE("Grid zero methods", "[grid]") {
    Grid g(2, 2, 2, 1.0, 1.0, 1.0);
    g.u[0] = 5.0; g.p[0] = 3.0;
    g.zero_velocity();
    for (double x : g.u) REQUIRE_THAT(x, WithinAbs(0.0, 1e-15));
    g.zero_pressure();
    for (double x : g.p) REQUIRE_THAT(x, WithinAbs(0.0, 1e-15));
}
