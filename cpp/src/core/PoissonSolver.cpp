#include "PoissonSolver.hpp"
#include <algorithm>
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
    // Non-uniform stencil: d²p/dx² at cell i uses:
    //   scale_r = 1/(dx_right * dx_cell),  dx_right = (dxv[i]+dxv[i+1])/2
    //   scale_l = 1/(dx_left  * dx_cell),  dx_left  = (dxv[i-1]+dxv[i])/2
    // For uniform grid (all dxv[i]=dx): scale = 1/dx², recovering original behavior.

    for (int i = 0; i < g.Nx; i++) {
        // Face-to-face distances (cell-center to cell-center)
        double hxl = (i > 0)      ? 0.5*(g.dxv[i-1]+g.dxv[i]) : g.dxv[i];
        double hxr = (i < g.Nx-1) ? 0.5*(g.dxv[i]+g.dxv[i+1]) : g.dxv[i];
        double sl_x = 1.0 / (hxl * g.dxv[i]);
        double sr_x = 1.0 / (hxr * g.dxv[i]);

        for (int j = 0; j < g.Ny; j++) {
            double hyl = (j > 0)      ? 0.5*(g.dyv[j-1]+g.dyv[j]) : g.dyv[j];
            double hyr = (j < g.Ny-1) ? 0.5*(g.dyv[j]+g.dyv[j+1]) : g.dyv[j];
            double sl_y = 1.0 / (hyl * g.dyv[j]);
            double sr_y = 1.0 / (hyr * g.dyv[j]);

            for (int k = 0; k < g.Nz; k++) {
                double hzl = (k > 0)      ? 0.5*(g.dzv[k-1]+g.dzv[k]) : g.dzv[k];
                double hzr = (k < g.Nz-1) ? 0.5*(g.dzv[k]+g.dzv[k+1]) : g.dzv[k];
                double sl_z = 1.0 / (hzl * g.dzv[k]);
                double sr_z = 1.0 / (hzr * g.dzv[k]);

                int c = g.c_idx(i,j,k);
                if (g.cell_type[c] == CellType::SOLID) {
                    Ax[c] = x[c];
                    continue;
                }

                double diag = 0.0;
                double off  = 0.0;

                // Lambda: process one neighbor face with its scale coefficient
                auto face = [&](int ni, int nj, int nk, double scale) {
                    if (ni < 0 || ni >= g.Nx ||
                        nj < 0 || nj >= g.Ny ||
                        nk < 0 || nk >= g.Nz) return;  // Neumann: zero contribution

                    int nc = g.c_idx(ni, nj, nk);
                    CellType ct = g.cell_type[nc];

                    if (ct == CellType::OUTFLOW) {
                        diag += 2.0 * scale;  // Dirichlet p=0
                    } else if (ct == CellType::SOLID) {
                        // Neumann ∂p/∂n=0: zero contribution
                    } else {
                        diag += scale;
                        off  += scale * x[nc];
                    }
                };

                face(i-1, j,   k,   sl_x);
                face(i+1, j,   k,   sr_x);
                face(i,   j-1, k,   sl_y);
                face(i,   j+1, k,   sr_y);
                face(i,   j,   k-1, sl_z);
                face(i,   j,   k+1, sr_z);

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

// ============================================================
//  MultigridPoissonSolver implementation
// ============================================================

MultigridPoissonSolver::MultigridPoissonSolver(PoissonSolverParams params)
    : params_(params) {}

// ---------------------------------------------------------------------------
// smooth: Gauss-Seidel point-by-point sweep.
// For each non-SOLID cell, update x[c] = (b[c] + sum_nb(scale_nb * x[nb])) / diag
// using the same non-uniform stencil as apply_laplacian().
// ---------------------------------------------------------------------------
void MultigridPoissonSolver::smooth(
        int Nx, int Ny, int Nz,
        const std::vector<double>& dxv,
        const std::vector<double>& dyv,
        const std::vector<double>& dzv,
        const std::vector<CellType>& cell_type,
        std::vector<double>& x,
        const std::vector<double>& b,
        int n_sweeps)
{
    auto cidx = [&](int i, int j, int k){ return i*(Ny*Nz) + j*Nz + k; };

    for (int sweep = 0; sweep < n_sweeps; ++sweep) {
        for (int i = 0; i < Nx; i++) {
            double hxl = (i > 0)      ? 0.5*(dxv[i-1]+dxv[i]) : dxv[i];
            double hxr = (i < Nx-1)   ? 0.5*(dxv[i]+dxv[i+1]) : dxv[i];
            double sl_x = 1.0 / (hxl * dxv[i]);
            double sr_x = 1.0 / (hxr * dxv[i]);

            for (int j = 0; j < Ny; j++) {
                double hyl = (j > 0)      ? 0.5*(dyv[j-1]+dyv[j]) : dyv[j];
                double hyr = (j < Ny-1)   ? 0.5*(dyv[j]+dyv[j+1]) : dyv[j];
                double sl_y = 1.0 / (hyl * dyv[j]);
                double sr_y = 1.0 / (hyr * dyv[j]);

                for (int k = 0; k < Nz; k++) {
                    int c = cidx(i,j,k);
                    if (cell_type[c] == CellType::SOLID) continue;

                    double hzl = (k > 0)      ? 0.5*(dzv[k-1]+dzv[k]) : dzv[k];
                    double hzr = (k < Nz-1)   ? 0.5*(dzv[k]+dzv[k+1]) : dzv[k];
                    double sl_z = 1.0 / (hzl * dzv[k]);
                    double sr_z = 1.0 / (hzr * dzv[k]);

                    double diag = 0.0;
                    double off  = 0.0;

                    auto face = [&](int ni, int nj, int nk, double scale) {
                        if (ni < 0 || ni >= Nx || nj < 0 || nj >= Ny ||
                            nk < 0 || nk >= Nz) return;
                        int nc = cidx(ni, nj, nk);
                        CellType ct = cell_type[nc];
                        if (ct == CellType::OUTFLOW) {
                            diag += 2.0 * scale;
                        } else if (ct == CellType::SOLID) {
                            // Neumann: zero contribution
                        } else {
                            diag += scale;
                            off  += scale * x[nc];
                        }
                    };

                    face(i-1, j,   k,   sl_x);
                    face(i+1, j,   k,   sr_x);
                    face(i,   j-1, k,   sl_y);
                    face(i,   j+1, k,   sr_y);
                    face(i,   j,   k-1, sl_z);
                    face(i,   j,   k+1, sr_z);

                    if (diag < 1e-30) diag = 1.0;
                    x[c] = (b[c] + off) / diag;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// compute_residual: r = b - A*x  (same stencil as apply_laplacian)
// ---------------------------------------------------------------------------
void MultigridPoissonSolver::compute_residual(
        int Nx, int Ny, int Nz,
        const std::vector<double>& dxv,
        const std::vector<double>& dyv,
        const std::vector<double>& dzv,
        const std::vector<CellType>& cell_type,
        const std::vector<double>& x,
        const std::vector<double>& b,
        std::vector<double>& r)
{
    auto cidx = [&](int i, int j, int k){ return i*(Ny*Nz) + j*Nz + k; };

    for (int i = 0; i < Nx; i++) {
        double hxl = (i > 0)    ? 0.5*(dxv[i-1]+dxv[i]) : dxv[i];
        double hxr = (i < Nx-1) ? 0.5*(dxv[i]+dxv[i+1]) : dxv[i];
        double sl_x = 1.0 / (hxl * dxv[i]);
        double sr_x = 1.0 / (hxr * dxv[i]);

        for (int j = 0; j < Ny; j++) {
            double hyl = (j > 0)    ? 0.5*(dyv[j-1]+dyv[j]) : dyv[j];
            double hyr = (j < Ny-1) ? 0.5*(dyv[j]+dyv[j+1]) : dyv[j];
            double sl_y = 1.0 / (hyl * dyv[j]);
            double sr_y = 1.0 / (hyr * dyv[j]);

            for (int k = 0; k < Nz; k++) {
                int c = cidx(i,j,k);
                if (cell_type[c] == CellType::SOLID) { r[c] = 0.0; continue; }

                double hzl = (k > 0)    ? 0.5*(dzv[k-1]+dzv[k]) : dzv[k];
                double hzr = (k < Nz-1) ? 0.5*(dzv[k]+dzv[k+1]) : dzv[k];
                double sl_z = 1.0 / (hzl * dzv[k]);
                double sr_z = 1.0 / (hzr * dzv[k]);

                double diag = 0.0;
                double off  = 0.0;

                auto face = [&](int ni, int nj, int nk, double scale) {
                    if (ni < 0 || ni >= Nx || nj < 0 || nj >= Ny ||
                        nk < 0 || nk >= Nz) return;
                    int nc = cidx(ni, nj, nk);
                    CellType ct = cell_type[nc];
                    if (ct == CellType::OUTFLOW) {
                        diag += 2.0 * scale;
                    } else if (ct == CellType::SOLID) {
                        // Neumann: zero
                    } else {
                        diag += scale;
                        off  += scale * x[nc];
                    }
                };

                face(i-1, j,   k,   sl_x);
                face(i+1, j,   k,   sr_x);
                face(i,   j-1, k,   sl_y);
                face(i,   j+1, k,   sr_y);
                face(i,   j,   k-1, sl_z);
                face(i,   j,   k+1, sr_z);

                if (diag < 1e-30) diag = 1.0;
                double Ax_c = diag * x[c] - off;
                r[c] = b[c] - Ax_c;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// restrict_residual: injection from fine (Nxf,Nyf,Nzf) to coarse (Nxc,Nyc,Nzc)
// ---------------------------------------------------------------------------
void MultigridPoissonSolver::restrict_residual(
        int Nxf, int Nyf, int Nzf,
        const std::vector<double>& r_fine,
        int Nxc, int Nyc, int Nzc,
        std::vector<double>& r_coarse)
{
    for (int ic = 0; ic < Nxc; ic++)
    for (int jc = 0; jc < Nyc; jc++)
    for (int kc = 0; kc < Nzc; kc++) {
        int if_ = 2*ic, jf = 2*jc, kf = 2*kc;
        if (if_ < Nxf && jf < Nyf && kf < Nzf)
            r_coarse[ic*Nyc*Nzc + jc*Nzc + kc] =
                r_fine[if_*Nyf*Nzf + jf*Nzf + kf];
        else
            r_coarse[ic*Nyc*Nzc + jc*Nzc + kc] = 0.0;
    }
}

// ---------------------------------------------------------------------------
// prolongate: nearest-neighbor injection from coarse to fine
// ---------------------------------------------------------------------------
void MultigridPoissonSolver::prolongate(
        int Nxc, int Nyc, int Nzc,
        const std::vector<double>& e_coarse,
        int Nxf, int Nyf, int Nzf,
        std::vector<double>& e_fine)
{
    std::fill(e_fine.begin(), e_fine.end(), 0.0);

    for (int ic = 0; ic < Nxc; ic++)
    for (int jc = 0; jc < Nyc; jc++)
    for (int kc = 0; kc < Nzc; kc++) {
        double ec = e_coarse[ic*Nyc*Nzc + jc*Nzc + kc];
        int if_  = std::min(2*ic,   Nxf-1);
        int if1  = std::min(2*ic+1, Nxf-1);
        int jf   = std::min(2*jc,   Nyf-1);
        int jf1  = std::min(2*jc+1, Nyf-1);
        int kf   = std::min(2*kc,   Nzf-1);
        int kf1  = std::min(2*kc+1, Nzf-1);
        // Apply to all 8 fine cells surrounding this coarse cell
        for (int ii : {if_, if1})
        for (int jj : {jf,  jf1})
        for (int kk : {kf,  kf1})
            e_fine[ii*Nyf*Nzf + jj*Nzf + kk] += ec;
    }
}

// ---------------------------------------------------------------------------
// v_cycle: recursive V-cycle
// ---------------------------------------------------------------------------
void MultigridPoissonSolver::v_cycle(
        int Nx, int Ny, int Nz,
        const std::vector<double>& dxv,
        const std::vector<double>& dyv,
        const std::vector<double>& dzv,
        const std::vector<CellType>& cell_type,
        std::vector<double>& x,
        const std::vector<double>& b,
        int level, int max_level)
{
    // Coarsest level: solve with many Gauss-Seidel sweeps
    if (level == max_level || (Nx <= 4 && Ny <= 4 && Nz <= 4)) {
        smooth(Nx, Ny, Nz, dxv, dyv, dzv, cell_type, x, b, 20);
        return;
    }

    // Pre-smoothing
    smooth(Nx, Ny, Nz, dxv, dyv, dzv, cell_type, x, b, 2);

    // Compute residual r = b - A*x
    int Nf = Nx * Ny * Nz;
    std::vector<double> res(Nf);
    compute_residual(Nx, Ny, Nz, dxv, dyv, dzv, cell_type, x, b, res);

    // Build coarse grid dimensions (halve each direction)
    int Nxc = std::max(1, Nx / 2);
    int Nyc = std::max(1, Ny / 2);
    int Nzc = std::max(1, Nz / 2);

    // Build coarse spacing arrays (combine pairs of fine cells)
    std::vector<double> dxvc(Nxc), dyvc(Nyc), dzvc(Nzc);
    for (int i = 0; i < Nxc; i++)
        dxvc[i] = dxv[2*i] + (2*i+1 < Nx ? dxv[2*i+1] : 0.0);
    for (int j = 0; j < Nyc; j++)
        dyvc[j] = dyv[2*j] + (2*j+1 < Ny ? dyv[2*j+1] : 0.0);
    for (int k = 0; k < Nzc; k++)
        dzvc[k] = dzv[2*k] + (2*k+1 < Nz ? dzv[2*k+1] : 0.0);

    // Coarse cell_type: use FLUID for all (approximate boundary handling)
    int Nc = Nxc * Nyc * Nzc;
    std::vector<CellType> cell_type_c(Nc, CellType::FLUID);

    // Restrict residual to coarse grid
    std::vector<double> r_coarse(Nc, 0.0);
    restrict_residual(Nx, Ny, Nz, res, Nxc, Nyc, Nzc, r_coarse);

    // Solve coarse error: e_c starts at zero
    std::vector<double> e_coarse(Nc, 0.0);
    v_cycle(Nxc, Nyc, Nzc, dxvc, dyvc, dzvc,
            cell_type_c, e_coarse, r_coarse,
            level + 1, max_level);

    // Prolongate error back to fine grid and correct x
    std::vector<double> e_fine(Nf, 0.0);
    prolongate(Nxc, Nyc, Nzc, e_coarse, Nx, Ny, Nz, e_fine);
    for (int idx = 0; idx < Nf; idx++) x[idx] += e_fine[idx];

    // Post-smoothing
    smooth(Nx, Ny, Nz, dxv, dyv, dzv, cell_type, x, b, 2);
}

// ---------------------------------------------------------------------------
// solve: run V-cycles until convergence or max_iter
// ---------------------------------------------------------------------------
double MultigridPoissonSolver::solve(Grid& grid, const std::vector<double>& rhs)
{
    int N = grid.Nx * grid.Ny * grid.Nz;
    std::vector<double>& x = grid.p;
    std::fill(x.begin(), x.end(), 0.0);

    // Compute rhs norm (ignore SOLID cells)
    double rhs_norm = 0.0;
    for (int i = 0; i < N; i++) {
        if (grid.cell_type[i] != CellType::SOLID)
            rhs_norm += rhs[i] * rhs[i];
    }
    rhs_norm = std::sqrt(rhs_norm);
    if (rhs_norm < 1e-15) {
        last_iter_ = 0; last_res_ = 0.0;
        return 0.0;
    }

    // Determine max multigrid level: stop when grid reaches ~4 cells per dim
    int min_dim = std::min({grid.Nx, grid.Ny, grid.Nz});
    int max_level = 0;
    {
        int d = min_dim;
        while (d > 4) { d /= 2; max_level++; }
        if (max_level < 1) max_level = 1;
    }

    std::vector<double> res_vec(N);

    for (int iter = 0; iter < params_.max_iter; iter++) {
        v_cycle(grid.Nx, grid.Ny, grid.Nz,
                grid.dxv, grid.dyv, grid.dzv,
                grid.cell_type, x, rhs,
                0, max_level);

        // Compute residual norm
        compute_residual(grid.Nx, grid.Ny, grid.Nz,
                         grid.dxv, grid.dyv, grid.dzv,
                         grid.cell_type, x, rhs, res_vec);
        double rr = 0.0;
        for (int i = 0; i < N; i++) rr += res_vec[i] * res_vec[i];
        double rel_res = std::sqrt(rr) / rhs_norm;

        if (rel_res < params_.tol) {
            last_iter_ = iter + 1;
            last_res_  = rel_res;
            return rel_res;
        }
    }

    // Recompute final residual
    compute_residual(grid.Nx, grid.Ny, grid.Nz,
                     grid.dxv, grid.dyv, grid.dzv,
                     grid.cell_type, x, rhs, res_vec);
    double rr = 0.0;
    for (int i = 0; i < N; i++) rr += res_vec[i] * res_vec[i];
    last_iter_ = params_.max_iter;
    last_res_  = std::sqrt(rr) / rhs_norm;
    return last_res_;
}

}  // namespace ffroom
