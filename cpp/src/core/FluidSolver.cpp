#include "FluidSolver.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace ffroom {

FluidSolver::FluidSolver(Grid& grid, BoundaryManager& bm, FluidSolverParams params)
    : grid_(grid), bm_(bm), params_(params),
      poisson_(params.poisson)
{
    u_prev_.resize(grid.u.size());
    v_prev_.resize(grid.v.size());
    w_prev_.resize(grid.w.size());

    // Initialize temperature field (T_initial everywhere; openings set to T_outside)
    std::fill(grid_.T.begin(), grid_.T.end(), params_.T_initial);
    if (params_.thermal) {
        bm_.apply_opening_temperature(grid_);
    }
}

// Interpolate u to cell center (i,j,k)
static inline double u_at_center(const Grid& g, int i, int j, int k) {
    return 0.5 * (g.u[g.u_idx(i,j,k)] + g.u[g.u_idx(i+1,j,k)]);
}
static inline double v_at_center(const Grid& g, int i, int j, int k) {
    return 0.5 * (g.v[g.v_idx(i,j,k)] + g.v[g.v_idx(i,j+1,k)]);
}
static inline double w_at_center(const Grid& g, int i, int j, int k) {
    return 0.5 * (g.w[g.w_idx(i,j,k)] + g.w[g.w_idx(i,j,k+1)]);
}

// First-order upwind advection of scalar phi transported by velocity uc
// Returns d(phi)/dt contribution from one direction
static inline double upwind1(double phi_m, double phi_c, double phi_p, double uc, double h) {
    if (uc > 0.0) return uc * (phi_c - phi_m) / h;
    else          return uc * (phi_p - phi_c) / h;
}

void FluidSolver::compute_intermediate() {
    const double dt  = params_.dt;
    const double nu  = params_.nu;
    const double idx = 1.0 / grid_.dx;
    const double idy = 1.0 / grid_.dy;
    const double idz = 1.0 / grid_.dz;
    const double idx2 = idx*idx, idy2 = idy*idy, idz2 = idz*idz;

    // --- u* for u-field (x-faces: i=0..Nx, j=0..Ny-1, k=0..Nz-1) ---
    for (int i = 1; i < grid_.Nx; i++)  // skip boundary faces i=0,Nx
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        int c = grid_.u_idx(i,j,k);
        double uc = grid_.u[c];

        // Advection: u transported by (u, v_at_uface, w_at_uface)
        // u at neighboring x-faces
        double u_im = (i > 1)        ? grid_.u[grid_.u_idx(i-1,j,k)] : 0.0;
        double u_ip = (i < grid_.Nx-1) ? grid_.u[grid_.u_idx(i+1,j,k)] : 0.0;

        // v interpolated to u-face position (i, j+0.5)
        double v_face = 0.0;
        if (j > 0 && j < grid_.Ny) {
            double v00 = grid_.v[grid_.v_idx(i-1, j, k)];
            double v10 = grid_.v[grid_.v_idx(i,   j, k)];
            v_face = 0.5 * (v00 + v10);
        }
        double u_jm = (j > 0)        ? grid_.u[grid_.u_idx(i,j-1,k)] : 0.0;
        double u_jp = (j < grid_.Ny-1) ? grid_.u[grid_.u_idx(i,j+1,k)] : 0.0;

        double w_face = 0.0;
        if (k > 0 && k < grid_.Nz) {
            double w00 = grid_.w[grid_.w_idx(i-1, j, k)];
            double w10 = grid_.w[grid_.w_idx(i,   j, k)];
            w_face = 0.5 * (w00 + w10);
        }
        double u_km = (k > 0)        ? grid_.u[grid_.u_idx(i,j,k-1)] : 0.0;
        double u_kp = (k < grid_.Nz-1) ? grid_.u[grid_.u_idx(i,j,k+1)] : 0.0;

        double adv = upwind1(u_im, uc, u_ip, uc,     grid_.dx)
                   + upwind1(u_jm, uc, u_jp, v_face, grid_.dy)
                   + upwind1(u_km, uc, u_kp, w_face, grid_.dz);

        // Diffusion: ν ∇²u
        double dif = nu * (idx2*(u_im - 2*uc + u_ip) +
                           idy2*(u_jm - 2*uc + u_jp) +
                           idz2*(u_km - 2*uc + u_kp));

        grid_.u_tmp[c] = uc + dt * (-adv + dif);
    }
    // Copy boundary faces unchanged
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        grid_.u_tmp[grid_.u_idx(0,      j,k)] = grid_.u[grid_.u_idx(0,      j,k)];
        grid_.u_tmp[grid_.u_idx(grid_.Nx,j,k)] = grid_.u[grid_.u_idx(grid_.Nx,j,k)];
    }

    // --- u* for v-field ---
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 1; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        int c = grid_.v_idx(i,j,k);
        double vc = grid_.v[c];

        double u_face = 0.0;
        if (i > 0 && i < grid_.Nx) {
            u_face = 0.5 * (grid_.u[grid_.u_idx(i, j-1, k)] + grid_.u[grid_.u_idx(i, j, k)]);
        }
        double v_im = (i > 0)        ? grid_.v[grid_.v_idx(i-1,j,k)] : 0.0;
        double v_ip = (i < grid_.Nx-1) ? grid_.v[grid_.v_idx(i+1,j,k)] : 0.0;
        double v_jm = (j > 1)        ? grid_.v[grid_.v_idx(i,j-1,k)] : 0.0;
        double v_jp = (j < grid_.Ny-1) ? grid_.v[grid_.v_idx(i,j+1,k)] : 0.0;

        double w_face = 0.0;
        if (k > 0 && k < grid_.Nz) {
            w_face = 0.5 * (grid_.w[grid_.w_idx(i, j-1, k)] + grid_.w[grid_.w_idx(i, j, k)]);
        }
        double v_km = (k > 0)        ? grid_.v[grid_.v_idx(i,j,k-1)] : 0.0;
        double v_kp = (k < grid_.Nz-1) ? grid_.v[grid_.v_idx(i,j,k+1)] : 0.0;

        double adv = upwind1(v_im, vc, v_ip, u_face, grid_.dx)
                   + upwind1(v_jm, vc, v_jp, vc,     grid_.dy)
                   + upwind1(v_km, vc, v_kp, w_face, grid_.dz);

        double dif = nu * (idx2*(v_im - 2*vc + v_ip) +
                           idy2*(v_jm - 2*vc + v_jp) +
                           idz2*(v_km - 2*vc + v_kp));

        grid_.v_tmp[c] = vc + dt * (-adv + dif);
    }
    for (int i = 0; i < grid_.Nx; i++)
    for (int k = 0; k < grid_.Nz; k++) {
        grid_.v_tmp[grid_.v_idx(i,0,      k)] = grid_.v[grid_.v_idx(i,0,      k)];
        grid_.v_tmp[grid_.v_idx(i,grid_.Ny,k)] = grid_.v[grid_.v_idx(i,grid_.Ny,k)];
    }

    // --- u* for w-field ---
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 1; k < grid_.Nz; k++) {
        int c = grid_.w_idx(i,j,k);
        double wc = grid_.w[c];

        double u_face = 0.0;
        if (i > 0 && i < grid_.Nx) {
            u_face = 0.5 * (grid_.u[grid_.u_idx(i, j, k-1)] + grid_.u[grid_.u_idx(i, j, k)]);
        }
        double w_im = (i > 0)        ? grid_.w[grid_.w_idx(i-1,j,k)] : 0.0;
        double w_ip = (i < grid_.Nx-1) ? grid_.w[grid_.w_idx(i+1,j,k)] : 0.0;

        double v_face = 0.0;
        if (j > 0 && j < grid_.Ny) {
            v_face = 0.5 * (grid_.v[grid_.v_idx(i, j, k-1)] + grid_.v[grid_.v_idx(i, j, k)]);
        }
        double w_jm = (j > 0)        ? grid_.w[grid_.w_idx(i,j-1,k)] : 0.0;
        double w_jp = (j < grid_.Ny-1) ? grid_.w[grid_.w_idx(i,j+1,k)] : 0.0;
        double w_km = (k > 1)        ? grid_.w[grid_.w_idx(i,j,k-1)] : 0.0;
        double w_kp = (k < grid_.Nz-1) ? grid_.w[grid_.w_idx(i,j,k+1)] : 0.0;

        double adv = upwind1(w_im, wc, w_ip, u_face, grid_.dx)
                   + upwind1(w_jm, wc, w_jp, v_face, grid_.dy)
                   + upwind1(w_km, wc, w_kp, wc,     grid_.dz);

        double dif = nu * (idx2*(w_im - 2*wc + w_ip) +
                           idy2*(w_jm - 2*wc + w_jp) +
                           idz2*(w_km - 2*wc + w_kp));

        grid_.w_tmp[c] = wc + dt * (-adv + dif);
    }
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++) {
        grid_.w_tmp[grid_.w_idx(i,j,0)]       = grid_.w[grid_.w_idx(i,j,0)];
        grid_.w_tmp[grid_.w_idx(i,j,grid_.Nz)] = grid_.w[grid_.w_idx(i,j,grid_.Nz)];
    }

    // Copy u* back to main arrays (BCs applied after)
    grid_.u = grid_.u_tmp;
    grid_.v = grid_.v_tmp;
    grid_.w = grid_.w_tmp;
}

void FluidSolver::build_poisson_rhs(std::vector<double>& rhs) const {
    // apply_laplacian computes L = -∇² (positive definite).
    // We solve: L*p = -(ρ/dt)*div(u*)  so that the standard correction
    // u = u* - (dt/ρ)*∇p reduces divergence.
    const double factor = -(params_.rho / params_.dt);
    const double idx = 1.0/grid_.dx, idy = 1.0/grid_.dy, idz = 1.0/grid_.dz;

    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        int c = grid_.c_idx(i,j,k);
        if (grid_.cell_type[c] == CellType::SOLID) {
            rhs[c] = 0.0;
            continue;
        }
        double div = (grid_.u[grid_.u_idx(i+1,j,k)] - grid_.u[grid_.u_idx(i,j,k)]) * idx
                   + (grid_.v[grid_.v_idx(i,j+1,k)] - grid_.v[grid_.v_idx(i,j,k)]) * idy
                   + (grid_.w[grid_.w_idx(i,j,k+1)] - grid_.w[grid_.w_idx(i,j,k)]) * idz;
        rhs[c] = factor * div;
    }
}

void FluidSolver::apply_correction() {
    const double dt_rho = params_.dt / params_.rho;
    const double idx = 1.0/grid_.dx, idy = 1.0/grid_.dy, idz = 1.0/grid_.dz;

    // Correct u: u -= (dt/rho) * dp/dx
    for (int i = 1; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        double dpdx = (grid_.p[grid_.c_idx(i,j,k)] - grid_.p[grid_.c_idx(i-1,j,k)]) * idx;
        grid_.u[grid_.u_idx(i,j,k)] -= dt_rho * dpdx;
    }

    // Correct v
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 1; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        double dpdy = (grid_.p[grid_.c_idx(i,j,k)] - grid_.p[grid_.c_idx(i,j-1,k)]) * idy;
        grid_.v[grid_.v_idx(i,j,k)] -= dt_rho * dpdy;
    }

    // Correct w
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 1; k < grid_.Nz; k++) {
        double dpdz = (grid_.p[grid_.c_idx(i,j,k)] - grid_.p[grid_.c_idx(i,j,k-1)]) * idz;
        grid_.w[grid_.w_idx(i,j,k)] -= dt_rho * dpdz;
    }
}

void FluidSolver::apply_buoyancy() {
    // Boussinesq: w* += g*β*(T_face - T_ref)*dt
    // T_face at z-face (i,j,k) = avg of cells k-1 and k.
    // Applies to interior z-faces k=1..Nz-1.
    const double coeff = params_.g_accel * params_.beta * params_.dt;
    const double T_ref = params_.T_initial;
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 1; k < grid_.Nz; k++) {
        double T_lo = grid_.T[grid_.c_idx(i, j, k-1)];
        double T_hi = grid_.T[grid_.c_idx(i, j, k)];
        grid_.w[grid_.w_idx(i, j, k)] += coeff * (0.5*(T_lo + T_hi) - T_ref);
    }
}

void FluidSolver::advect_temperature() {
    const double dt    = params_.dt;
    const double idx2  = 1.0/(grid_.dx*grid_.dx);
    const double idy2  = 1.0/(grid_.dy*grid_.dy);
    const double idz2  = 1.0/(grid_.dz*grid_.dz);
    const double alpha = params_.k_thermal / (params_.rho * params_.cp);

    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        int c = grid_.c_idx(i,j,k);
        if (grid_.cell_type[c] != CellType::FLUID) {
            grid_.T_tmp[c] = grid_.T[c];
            continue;
        }

        double Tc = grid_.T[c];
        double uc = u_at_center(grid_, i, j, k);
        double vc = v_at_center(grid_, i, j, k);
        double wc = w_at_center(grid_, i, j, k);

        // Neighbor T: Neumann (T=Tc) at domain boundary and SOLID walls
        auto get_T = [&](int ni, int nj, int nk) -> double {
            if (ni < 0 || ni >= grid_.Nx || nj < 0 || nj >= grid_.Ny ||
                nk < 0 || nk >= grid_.Nz) return Tc;
            int nc = grid_.c_idx(ni, nj, nk);
            return (grid_.cell_type[nc] == CellType::SOLID) ? Tc : grid_.T[nc];
        };

        double T_xm = get_T(i-1,j,k), T_xp = get_T(i+1,j,k);
        double T_ym = get_T(i,j-1,k), T_yp = get_T(i,j+1,k);
        double T_zm = get_T(i,j,k-1), T_zp = get_T(i,j,k+1);

        double adv = upwind1(T_xm, Tc, T_xp, uc, grid_.dx)
                   + upwind1(T_ym, Tc, T_yp, vc, grid_.dy)
                   + upwind1(T_zm, Tc, T_zp, wc, grid_.dz);

        double dif = alpha * (idx2*(T_xm - 2*Tc + T_xp) +
                              idy2*(T_ym - 2*Tc + T_yp) +
                              idz2*(T_zm - 2*Tc + T_zp));

        grid_.T_tmp[c] = Tc + dt * (-adv + dif);
    }
    std::swap(grid_.T, grid_.T_tmp);
}

double FluidSolver::compute_T_mean() const {
    double sum = 0.0;
    int    n   = 0;
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        int c = grid_.c_idx(i,j,k);
        if (grid_.cell_type[c] == CellType::FLUID) { sum += grid_.T[c]; ++n; }
    }
    return n > 0 ? sum / n : params_.T_initial;
}

double FluidSolver::compute_divergence_max() const {
    const double idx = 1.0/grid_.dx, idy = 1.0/grid_.dy, idz = 1.0/grid_.dz;
    double max_div = 0.0;
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        CellType ct = grid_.cell_type[grid_.c_idx(i,j,k)];
        // Only check pure FLUID cells; boundary cells have prescribed divergence
        if (ct != CellType::FLUID) continue;
        double div = std::abs(
            (grid_.u[grid_.u_idx(i+1,j,k)] - grid_.u[grid_.u_idx(i,j,k)]) * idx
          + (grid_.v[grid_.v_idx(i,j+1,k)] - grid_.v[grid_.v_idx(i,j,k)]) * idy
          + (grid_.w[grid_.w_idx(i,j,k+1)] - grid_.w[grid_.w_idx(i,j,k)]) * idz);
        max_div = std::max(max_div, div);
    }
    return max_div;
}

StepResult FluidSolver::step() {
    u_prev_ = grid_.u;
    v_prev_ = grid_.v;
    w_prev_ = grid_.w;

    // 1. Intermediate velocity (advection + diffusion)
    compute_intermediate();

    // 2. Buoyancy force on w* (must precede Poisson to project correctly)
    if (params_.buoyancy) apply_buoyancy();

    // 3. Enforce BCs on u*
    bm_.apply_noslip(grid_);
    bm_.apply_inflow(grid_);

    // 4. Solve pressure Poisson
    int N = grid_.Nx * grid_.Ny * grid_.Nz;
    std::vector<double> rhs(N);
    build_poisson_rhs(rhs);
    double p_res = poisson_.solve(grid_, rhs);

    // 5. Correct velocity
    apply_correction();

    // 6. Re-enforce BCs on corrected velocity
    bm_.apply_noslip(grid_);
    bm_.apply_inflow(grid_);
    bm_.apply_outflow(grid_);

    // 7. Temperature update (uses divergence-free corrected velocity)
    double T_mean = params_.T_initial;
    if (params_.thermal) {
        bm_.apply_opening_temperature(grid_);
        advect_temperature();
        bm_.apply_opening_temperature(grid_);
        T_mean = compute_T_mean();
    }

    // Diagnostics
    double vel_change = 0.0;
    for (size_t i = 0; i < grid_.u.size(); i++)
        vel_change = std::max(vel_change, std::abs(grid_.u[i] - u_prev_[i]));
    for (size_t i = 0; i < grid_.v.size(); i++)
        vel_change = std::max(vel_change, std::abs(grid_.v[i] - v_prev_[i]));
    for (size_t i = 0; i < grid_.w.size(); i++)
        vel_change = std::max(vel_change, std::abs(grid_.w[i] - w_prev_[i]));
    vel_change /= params_.dt;

    step_++;

    // Convergence: temperature-target or velocity steady-state
    bool converged;
    if (params_.thermal && params_.T_target < 1e29) {
        // Stop when mean room temperature reaches the target
        converged = (params_.T_target >= params_.T_initial)
                  ? (T_mean >= params_.T_target)
                  : (T_mean <= params_.T_target);
    } else {
        converged = (vel_change < params_.convergence_tol);
    }

    return {step_, vel_change, p_res, compute_divergence_max(), converged, T_mean};
}

StepResult FluidSolver::run(std::function<void(const StepResult&)> callback) {
    StepResult result{};
    for (int s = 0; s < params_.max_steps; s++) {
        result = step();
        if (callback) callback(result);
        if (result.converged) break;
    }
    return result;
}

}  // namespace ffroom
