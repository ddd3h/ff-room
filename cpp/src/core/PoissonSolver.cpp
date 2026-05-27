#include "PoissonSolver.hpp"
#include <cmath>
#include <stdexcept>

namespace ffroom {

PoissonSolver::PoissonSolver(PoissonSolverParams params)
    : params_(params) {}

// Computes Ax = L*x where L = -∇² (positive (semi)definite).
//
// [L*x]_c = sum_faces (x_c - x_nb) / h²
//
// Boundary treatment per face:
//   Domain wall / SOLID  → Neumann ∂p/∂n=0: ghost = x_c  → contribution = 0
//   OUTFLOW               → Dirichlet p=0:   ghost = -x_c → contribution = 2/h² * x_c
//   INFLOW / FLUID        → interior stencil: 1/h² * (x_c - x_nb)
void PoissonSolver::apply_laplacian(const Grid& g, const std::vector<double>& x,
                                     std::vector<double>& Ax) {
    const double idx2 = 1.0 / (g.dx * g.dx);
    const double idy2 = 1.0 / (g.dy * g.dy);
    const double idz2 = 1.0 / (g.dz * g.dz);

    for (int i = 0; i < g.Nx; i++) {
        for (int j = 0; j < g.Ny; j++) {
            for (int k = 0; k < g.Nz; k++) {
                int c = g.c_idx(i,j,k);
                if (g.cell_type[c] == CellType::SOLID) {
                    Ax[c] = x[c];  // identity row; no Poisson equation in solid
                    continue;
                }

                double diag = 0.0;
                double off  = 0.0;

                auto face = [&](int ni, int nj, int nk, double scale) {
                    // Out-of-domain → Neumann (contribution = 0)
                    if (ni < 0 || ni >= g.Nx ||
                        nj < 0 || nj >= g.Ny ||
                        nk < 0 || nk >= g.Nz) return;

                    int nc = g.c_idx(ni, nj, nk);
                    CellType ct = g.cell_type[nc];

                    if (ct == CellType::OUTFLOW) {
                        // Dirichlet p=0 at face: ghost = -x[c] → contribution = 2/h²
                        diag += 2.0 * scale;
                    } else if (ct == CellType::SOLID) {
                        // Neumann ∂p/∂n=0: ghost = x[c] → contribution = 0
                    } else {
                        // FLUID or INFLOW: standard interior stencil
                        diag += scale;
                        off  += scale * x[nc];
                    }
                };

                face(i-1, j,   k,   idx2);
                face(i+1, j,   k,   idx2);
                face(i,   j-1, k,   idy2);
                face(i,   j+1, k,   idy2);
                face(i,   j,   k-1, idz2);
                face(i,   j,   k+1, idz2);

                // Guard: if all faces are Neumann (isolated cell), make diagonal = 1
                // to avoid singular row in CG.
                if (diag < 1e-30) diag = 1.0;

                Ax[c] = diag * x[c] - off;
            }
        }
    }
}

double PoissonSolver::solve(Grid& grid, const std::vector<double>& rhs) {
    int N = grid.Nx * grid.Ny * grid.Nz;
    std::vector<double>& x = grid.p;
    std::fill(x.begin(), x.end(), 0.0);  // initial guess = 0

    std::vector<double> Ax(N, 0.0);
    std::vector<double> r(N), d(N), Ad(N);

    // r = rhs - A*x  (x=0 → r = rhs initially; zero out solid/outflow)
    for (int i = 0; i < N; i++) {
        CellType ct = grid.cell_type[i];
        r[i] = (ct == CellType::SOLID) ? 0.0 : rhs[i];
        d[i] = r[i];
    }

    double rr = 0.0;
    for (int i = 0; i < N; i++) rr += r[i]*r[i];
    double rhs_norm = std::sqrt(rr);
    if (rhs_norm < 1e-15) {
        last_iter_ = 0; last_res_ = 0.0;
        return 0.0;
    }

    for (int iter = 0; iter < params_.max_iter; iter++) {
        apply_laplacian(grid, d, Ad);

        double dAd = 0.0;
        for (int i = 0; i < N; i++) dAd += d[i]*Ad[i];
        if (std::abs(dAd) < 1e-30) break;

        double alpha = rr / dAd;
        double rr_new = 0.0;
        for (int i = 0; i < N; i++) {
            x[i] += alpha * d[i];
            r[i] -= alpha * Ad[i];
            rr_new += r[i]*r[i];
        }

        double res = std::sqrt(rr_new) / rhs_norm;
        if (res < params_.tol) {
            last_iter_ = iter + 1;
            last_res_  = res;
            return res;
        }

        double beta = rr_new / rr;
        for (int i = 0; i < N; i++) d[i] = r[i] + beta * d[i];
        rr = rr_new;
    }

    last_iter_ = params_.max_iter;
    last_res_  = std::sqrt(rr) / rhs_norm;
    return last_res_;
}

}  // namespace ffroom
