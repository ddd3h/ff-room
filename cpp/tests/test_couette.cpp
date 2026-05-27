#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/Grid.hpp"
#include "core/FluidSolver.hpp"
#include "core/BoundaryManager.hpp"
#include <cmath>

using namespace ffroom;
using Catch::Matchers::WithinAbs;

// Couette flow benchmark:
// 2D-like flow in z-direction: u = U_wall * z/H, v=0, w=0, p=const
// Driven by moving top wall (z=H), static bottom (z=0).
// Steady-state analytic: u(z) = U_wall * z / H
//
// Implementation: thin domain in y (Ny=2), drive by inflow at top face.
// Steady state should match linear profile to within ~5% on coarse grid.

TEST_CASE("Couette flow: linear profile converged", "[couette][benchmark]") {
    const int   Nx = 4, Ny = 4, Nz = 10;
    const double Lx=1.0, Ly=1.0, Lz=1.0;
    const double U_wall = 1.0;  // [m/s] top-wall x-velocity

    Grid g(Nx, Ny, Nz, Lx, Ly, Lz);
    BoundaryManager bm;

    // Drive: set inflow on top x-faces of all cells at z=Nz-1 (ceiling)
    // Instead, we directly impose boundary via fan (set u on ceiling to U_wall)
    // Simple approach: manually set top u-faces and run.
    // Actually: Couette driven by moving wall. In no-slip BM, ceiling velocity=U_wall.
    // We override apply_noslip by adding a "moving wall" fan at the ceiling.
    FanBC top_wall;
    top_wall.i_min = 0; top_wall.i_max = Nx-1;
    top_wall.j_min = 0; top_wall.j_max = Ny-1;
    top_wall.k_min = Nz-1; top_wall.k_max = Nz-1;
    top_wall.vel   = {U_wall, 0.0, 0.0};
    bm.add_fan(g, top_wall);

    FluidSolverParams params;
    params.dt             = 0.005;
    params.max_steps      = 2000;
    params.convergence_tol = 1e-5;
    params.nu             = 1.5e-5;
    params.poisson.tol    = 1e-7;

    FluidSolver solver(g, bm, params);
    auto result = solver.run();

    // Extract u profile at x=Nx/2, y=Ny/2 (cell centers interpolated from u-faces)
    int i = Nx/2, j = Ny/2;
    double max_err = 0.0;
    for (int k = 0; k < Nz; k++) {
        double z_center = (k + 0.5) * g.dz;  // cell center z
        double u_sim  = 0.5 * (g.u[g.u_idx(i,j,k)] + g.u[g.u_idx(i+1,j,k)]);
        double u_anal = U_wall * z_center / Lz;
        max_err = std::max(max_err, std::abs(u_sim - u_anal));
    }

    // Coarse grid (Nz=10): expect < 15% relative error
    REQUIRE(max_err < 0.15 * U_wall);

    // Also check convergence happened
    REQUIRE(result.converged);
}
