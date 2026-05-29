#pragma once
#include "Grid.hpp"

namespace ffroom {

// Solves: ∇²p = rhs
// where rhs = (rho/dt) * divergence(u*)
// Boundary: dp/dn = 0 on no-slip walls (Neumann)
//           p = 0 on outflow faces (Dirichlet)
//
// Method: Conjugate Gradient (CG) with Jacobi preconditioner
// Convergence: ||r||_2 / ||rhs||_2 < tol

struct PoissonSolverParams {
    int    max_iter   = 1000;
    double tol        = 1e-6;
    bool   use_openmp = false;
};

class PoissonSolver {
public:
    explicit PoissonSolver(PoissonSolverParams params = {});

    // Fills grid.p. Returns residual at convergence.
    double solve(Grid& grid, const std::vector<double>& rhs);

    int    last_iterations() const { return last_iter_; }
    double last_residual()   const { return last_res_; }

private:
    PoissonSolverParams params_;
    int    last_iter_ = 0;
    double last_res_  = 0.0;

    // Apply Laplacian stencil to p, result stored in Ap
    void apply_laplacian(const Grid& grid, const std::vector<double>& x,
                         std::vector<double>& Ax);
};

// Multigrid V-cycle solver for the pressure Poisson equation.
// O(N) convergence as an alternative to CG (O(N^1.5)).
class MultigridPoissonSolver {
public:
    explicit MultigridPoissonSolver(PoissonSolverParams params = {});

    // Same interface as PoissonSolver — fills grid.p, returns residual.
    double solve(Grid& grid, const std::vector<double>& rhs);

    int    last_iterations() const { return last_iter_; }
    double last_residual()   const { return last_res_; }

private:
    PoissonSolverParams params_;
    int    last_iter_ = 0;
    double last_res_  = 0.0;

    // V-cycle: one multigrid iteration (recursive coarse-grid correction)
    void v_cycle(int Nx, int Ny, int Nz,
                 const std::vector<double>& dxv,
                 const std::vector<double>& dyv,
                 const std::vector<double>& dzv,
                 const std::vector<CellType>& cell_type,
                 std::vector<double>& x,
                 const std::vector<double>& b,
                 int level, int max_level);

    // Gauss-Seidel smoother (simple sweep)
    void smooth(int Nx, int Ny, int Nz,
                const std::vector<double>& dxv,
                const std::vector<double>& dyv,
                const std::vector<double>& dzv,
                const std::vector<CellType>& cell_type,
                std::vector<double>& x,
                const std::vector<double>& b,
                int n_sweeps);

    // Restriction: fine -> coarse (injection at every other point)
    void restrict_residual(int Nxf, int Nyf, int Nzf,
                           const std::vector<double>& r_fine,
                           int Nxc, int Nyc, int Nzc,
                           std::vector<double>& r_coarse);

    // Prolongation: coarse -> fine (nearest-neighbor trilinear)
    void prolongate(int Nxc, int Nyc, int Nzc,
                    const std::vector<double>& e_coarse,
                    int Nxf, int Nyf, int Nzf,
                    std::vector<double>& e_fine);

    // Compute residual r = b - A*x
    void compute_residual(int Nx, int Ny, int Nz,
                          const std::vector<double>& dxv,
                          const std::vector<double>& dyv,
                          const std::vector<double>& dzv,
                          const std::vector<CellType>& cell_type,
                          const std::vector<double>& x,
                          const std::vector<double>& b,
                          std::vector<double>& r);
};

}  // namespace ffroom
