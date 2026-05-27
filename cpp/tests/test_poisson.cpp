#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/Grid.hpp"
#include "core/PoissonSolver.hpp"
#include <cmath>
#include <set>

using namespace ffroom;
using Catch::Matchers::WithinAbs;

// Verify CG solves a known RHS on small grid (all-Neumann case).
// Manufactured solution: p = cos(πx/L) * cos(πy/L) * cos(πz/L)
// Lap(p) = -3*(π/L)^2 * p
// With Neumann BCs dp/dn=0 on all walls -> pure Neumann, singular system.
// We pin one cell to zero to get unique solution.

TEST_CASE("Poisson CG converges on small grid", "[poisson]") {
    int N = 8;
    Grid g(N, N, N, 1.0, 1.0, 1.0);

    // All cells OUTFLOW so p=0 on all boundaries -> Dirichlet system
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        bool boundary = (i==0||i==N-1||j==0||j==N-1||k==0||k==N-1);
        if (boundary) g.cell_type[g.c_idx(i,j,k)] = CellType::OUTFLOW;
    }

    // Build RHS from known p = sin(πx)*sin(πy)*sin(πz) on unit cube interior
    // Lap(p) = -3π² p
    double dx = g.dx, dy = g.dy, dz = g.dz;
    std::vector<double> rhs(N*N*N, 0.0);
    std::vector<double> p_exact(N*N*N, 0.0);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        double x = (i+0.5)*dx, y = (j+0.5)*dy, z = (k+0.5)*dz;
        double pe = std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
        p_exact[g.c_idx(i,j,k)] = pe;
        bool boundary = (i==0||i==N-1||j==0||j==N-1||k==0||k==N-1);
        if (!boundary)
            rhs[g.c_idx(i,j,k)] = -3.0*M_PI*M_PI * pe;
    }

    PoissonSolverParams params;
    params.tol = 1e-8;
    PoissonSolver solver(params);
    double res = solver.solve(g, rhs);

    REQUIRE(res < params.tol * 10);  // allow slight overshoot

    // Check solution accuracy in interior
    double max_err = 0.0;
    for (int i = 1; i < N-1; i++)
    for (int j = 1; j < N-1; j++)
    for (int k = 1; k < N-1; k++) {
        int c = g.c_idx(i,j,k);
        max_err = std::max(max_err, std::abs(g.p[c] - p_exact[c]));
    }
    // First-order FD on 8-cell grid: expect ~1% error
    REQUIRE(max_err < 0.05);
}
