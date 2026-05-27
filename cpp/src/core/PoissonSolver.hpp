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
    int    max_iter = 1000;
    double tol      = 1e-6;
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

}  // namespace ffroom
