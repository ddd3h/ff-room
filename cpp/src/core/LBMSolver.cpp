#include "LBMSolver.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace ffroom {

// ---------------------------------------------------------------------------
// D3Q19 lattice constants
// ---------------------------------------------------------------------------

// Velocity vectors for each of the 19 directions
const int LBMSolver::EX[19] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0};
const int LBMSolver::EY[19] = { 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1};
const int LBMSolver::EZ[19] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1};

// Weights
const double LBMSolver::W[19] = {
    1.0/3.0,
    1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
};

// Opposite direction: OPP[i] such that E[OPP[i]] = -E[i]
const int LBMSolver::OPP[19] = {0, 2,1, 4,3, 6,5, 8,7, 10,9, 12,11, 14,13, 16,15, 18,17};

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

LBMSolver::LBMSolver(Grid& grid, BoundaryManager& bm, double tau,
                     int max_steps, double convergence_tol, double dt)
    : grid_(grid), bm_(bm), tau_(tau),
      max_steps_(max_steps), convergence_tol_(convergence_tol), dt_(dt)
{
    // Use the mean grid spacing as the reference length (dx).
    // lattice velocity = physical velocity * dt / dx_ref
    // A uniform dx across all directions keeps things simple.
    // We use grid_.dx (= Lx/Nx) as the reference dx.
    double dx_ref = grid_.dx;
    phys_to_lat_ = dt_ / dx_ref;         // [s/m] * [m] = dimensionless
    lat_to_phys_ = dx_ref / dt_;         // [m/s]

    int N = grid_.Nx * grid_.Ny * grid_.Nz;
    f_.resize(Q * N, 0.0);
    f_tmp_.resize(Q * N, 0.0);
    ux_.resize(N, 0.0);
    uy_.resize(N, 0.0);
    uz_.resize(N, 0.0);
    rho_.resize(N, 1.0);

    u_prev_.resize(grid_.u.size(), 0.0);
    v_prev_.resize(grid_.v.size(), 0.0);
    w_prev_.resize(grid_.w.size(), 0.0);

    initialize_f();
}

// ---------------------------------------------------------------------------
// Initialize distributions to equilibrium at rest (rho=1, u=0)
// ---------------------------------------------------------------------------

void LBMSolver::initialize_f() {
    int N = grid_.Nx * grid_.Ny * grid_.Nz;
    for (int q = 0; q < Q; q++) {
        for (int n = 0; n < N; n++) {
            f_[q * N + n] = W[q];  // feq(rho=1, u=0) = W[q]
        }
    }
    // Apply inflow BC to initialize inflow cells
    apply_inflow_bc();
}

// ---------------------------------------------------------------------------
// Inflow BC: set distributions to equilibrium with the fan velocity
// ---------------------------------------------------------------------------

void LBMSolver::apply_inflow_bc() {
    const auto& fans = bm_.fans();
    for (const auto& fan : fans) {
        // Convert physical velocity [m/s] to lattice velocity [dimensionless]
        double vx = fan.vel[0] * phys_to_lat_;
        double vy = fan.vel[1] * phys_to_lat_;
        double vz = fan.vel[2] * phys_to_lat_;
        for (int i = fan.i_min; i <= fan.i_max; i++)
        for (int j = fan.j_min; j <= fan.j_max; j++)
        for (int k = fan.k_min; k <= fan.k_max; k++) {
            if (i < 0 || i >= grid_.Nx || j < 0 || j >= grid_.Ny || k < 0 || k >= grid_.Nz)
                continue;
            int c = grid_.c_idx(i, j, k);
            if (grid_.cell_type[c] != CellType::INFLOW) continue;
            for (int q = 0; q < Q; q++) {
                f_[f_idx(q, i, j, k)] = feq(q, 1.0, vx, vy, vz);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// BGK collision: f_i += -(f_i - f_i^eq) / tau
// Skip SOLID cells (they only participate via bounce-back during streaming)
// ---------------------------------------------------------------------------

void LBMSolver::collision() {
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        int c = grid_.c_idx(i, j, k);
        CellType ct = grid_.cell_type[c];
        if (ct == CellType::SOLID) continue;  // walls: pure bounce-back, no collision

        // Compute macroscopic quantities from current f
        double rho = 0.0, ux = 0.0, uy = 0.0, uz = 0.0;
        for (int q = 0; q < Q; q++) {
            double fq = f_[f_idx(q, i, j, k)];
            rho += fq;
            ux  += fq * EX[q];
            uy  += fq * EY[q];
            uz  += fq * EZ[q];
        }
        if (rho > 1e-14) {
            ux /= rho;
            uy /= rho;
            uz /= rho;
        }

        // For INFLOW cells, force equilibrium at prescribed velocity
        if (ct == CellType::INFLOW) {
            // Already set to eq in apply_inflow_bc(); skip collision here to
            // keep the distribution exactly at eq (prevents drift)
            continue;
        }

        // BGK collision
        double inv_tau = 1.0 / tau_;
        for (int q = 0; q < Q; q++) {
            int idx = f_idx(q, i, j, k);
            double f_eq = feq(q, rho, ux, uy, uz);
            f_[idx] += inv_tau * (f_eq - f_[idx]);
        }
    }
}

// ---------------------------------------------------------------------------
// Streaming: push post-collision distributions to neighbors.
// Bounce-back for solid walls and domain boundaries.
// ---------------------------------------------------------------------------

void LBMSolver::streaming() {
    // Initialize f_tmp with zeros (unfilled entries stay zero)
    std::fill(f_tmp_.begin(), f_tmp_.end(), 0.0);

    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        int c = grid_.c_idx(i, j, k);
        CellType ct = grid_.cell_type[c];
        if (ct == CellType::SOLID) continue;  // SOLID cells don't stream out

        for (int q = 0; q < Q; q++) {
            double fq = f_[f_idx(q, i, j, k)];

            int ni = i + EX[q];
            int nj = j + EY[q];
            int nk = k + EZ[q];

            // Check if neighbor is in-domain and not SOLID
            bool out_of_bounds = (ni < 0 || ni >= grid_.Nx ||
                                  nj < 0 || nj >= grid_.Ny ||
                                  nk < 0 || nk >= grid_.Nz);
            bool is_solid_neighbor = false;
            if (!out_of_bounds) {
                is_solid_neighbor = (grid_.cell_type[grid_.c_idx(ni, nj, nk)] == CellType::SOLID);
            }

            if (out_of_bounds || is_solid_neighbor) {
                // Bounce-back: reflect back to current cell in opposite direction
                f_tmp_[f_idx(OPP[q], i, j, k)] += fq;
            } else {
                // Normal stream
                f_tmp_[f_idx(q, ni, nj, nk)] += fq;
            }
        }
    }

    // Swap buffers
    std::swap(f_, f_tmp_);
}

// ---------------------------------------------------------------------------
// Outflow BC: zero-gradient — copy f from nearest interior neighbor
// ---------------------------------------------------------------------------

void LBMSolver::apply_outflow_bc() {
    const auto& outflows = bm_.outflows();
    for (const auto& obc : outflows) {
        if (obc.axis == 0) {
            // x-axis outflow
            int i_oflow = (obc.side == 0) ? 0 : grid_.Nx - 1;
            int i_src   = (obc.side == 0) ? 1 : grid_.Nx - 2;
            if (i_src < 0 || i_src >= grid_.Nx) continue;
            for (int j = obc.a_min; j <= obc.a_max; j++)
            for (int k = obc.b_min; k <= obc.b_max; k++) {
                if (j < 0 || j >= grid_.Ny || k < 0 || k >= grid_.Nz) continue;
                int co = grid_.c_idx(i_oflow, j, k);
                if (grid_.cell_type[co] != CellType::OUTFLOW) continue;
                for (int q = 0; q < Q; q++)
                    f_[f_idx(q, i_oflow, j, k)] = f_[f_idx(q, i_src, j, k)];
            }
        } else if (obc.axis == 1) {
            // y-axis outflow
            int j_oflow = (obc.side == 0) ? 0 : grid_.Ny - 1;
            int j_src   = (obc.side == 0) ? 1 : grid_.Ny - 2;
            if (j_src < 0 || j_src >= grid_.Ny) continue;
            for (int i = obc.a_min; i <= obc.a_max; i++)
            for (int k = obc.b_min; k <= obc.b_max; k++) {
                if (i < 0 || i >= grid_.Nx || k < 0 || k >= grid_.Nz) continue;
                int co = grid_.c_idx(i, j_oflow, k);
                if (grid_.cell_type[co] != CellType::OUTFLOW) continue;
                for (int q = 0; q < Q; q++)
                    f_[f_idx(q, i, j_oflow, k)] = f_[f_idx(q, i, j_src, k)];
            }
        } else {
            // z-axis outflow
            int k_oflow = (obc.side == 0) ? 0 : grid_.Nz - 1;
            int k_src   = (obc.side == 0) ? 1 : grid_.Nz - 2;
            if (k_src < 0 || k_src >= grid_.Nz) continue;
            for (int i = obc.a_min; i <= obc.a_max; i++)
            for (int j = obc.b_min; j <= obc.b_max; j++) {
                if (i < 0 || i >= grid_.Nx || j < 0 || j >= grid_.Ny) continue;
                int co = grid_.c_idx(i, j, k_oflow);
                if (grid_.cell_type[co] != CellType::OUTFLOW) continue;
                for (int q = 0; q < Q; q++)
                    f_[f_idx(q, i, j, k_oflow)] = f_[f_idx(q, i, j, k_src)];
            }
        }
    }

    // Also handle openings (same structure as OutflowBC)
    const auto& openings = bm_.openings();
    for (const auto& obc : openings) {
        if (obc.axis == 0) {
            int i_oflow = (obc.side == 0) ? 0 : grid_.Nx - 1;
            int i_src   = (obc.side == 0) ? 1 : grid_.Nx - 2;
            if (i_src < 0 || i_src >= grid_.Nx) continue;
            for (int j = obc.a_min; j <= obc.a_max; j++)
            for (int k = obc.b_min; k <= obc.b_max; k++) {
                if (j < 0 || j >= grid_.Ny || k < 0 || k >= grid_.Nz) continue;
                int co = grid_.c_idx(i_oflow, j, k);
                if (grid_.cell_type[co] != CellType::OUTFLOW) continue;
                for (int q = 0; q < Q; q++)
                    f_[f_idx(q, i_oflow, j, k)] = f_[f_idx(q, i_src, j, k)];
            }
        } else if (obc.axis == 1) {
            int j_oflow = (obc.side == 0) ? 0 : grid_.Ny - 1;
            int j_src   = (obc.side == 0) ? 1 : grid_.Ny - 2;
            if (j_src < 0 || j_src >= grid_.Ny) continue;
            for (int i = obc.a_min; i <= obc.a_max; i++)
            for (int k = obc.b_min; k <= obc.b_max; k++) {
                if (i < 0 || i >= grid_.Nx || k < 0 || k >= grid_.Nz) continue;
                int co = grid_.c_idx(i, j_oflow, k);
                if (grid_.cell_type[co] != CellType::OUTFLOW) continue;
                for (int q = 0; q < Q; q++)
                    f_[f_idx(q, i, j_oflow, k)] = f_[f_idx(q, i, j_src, k)];
            }
        } else {
            int k_oflow = (obc.side == 0) ? 0 : grid_.Nz - 1;
            int k_src   = (obc.side == 0) ? 1 : grid_.Nz - 2;
            if (k_src < 0 || k_src >= grid_.Nz) continue;
            for (int i = obc.a_min; i <= obc.a_max; i++)
            for (int j = obc.b_min; j <= obc.b_max; j++) {
                if (i < 0 || i >= grid_.Nx || j < 0 || j >= grid_.Ny) continue;
                int co = grid_.c_idx(i, j, k_oflow);
                if (grid_.cell_type[co] != CellType::OUTFLOW) continue;
                for (int q = 0; q < Q; q++)
                    f_[f_idx(q, i, j, k_oflow)] = f_[f_idx(q, i, j, k_src)];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Extract macroscopic fields from f; write velocities/pressure to grid
// ---------------------------------------------------------------------------

void LBMSolver::extract_macroscopic() {
    const int Nx = grid_.Nx;
    const int Ny = grid_.Ny;
    const int Nz = grid_.Nz;

    // Compute cell-centered rho, ux, uy, uz from distributions
    // Velocities are in lattice units; convert to physical [m/s] via lat_to_phys_
    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++) {
        int c = c_idx(i, j, k);
        CellType ct = grid_.cell_type[c];
        if (ct == CellType::SOLID) {
            rho_[c] = 1.0;
            ux_[c] = uy_[c] = uz_[c] = 0.0;
            continue;
        }

        double rho = 0.0, ux = 0.0, uy = 0.0, uz = 0.0;
        for (int q = 0; q < Q; q++) {
            double fq = f_[f_idx(q, i, j, k)];
            rho += fq;
            ux  += fq * EX[q];
            uy  += fq * EY[q];
            uz  += fq * EZ[q];
        }
        if (rho > 1e-14) {
            ux /= rho;
            uy /= rho;
            uz /= rho;
        }
        rho_[c] = rho;
        // Convert lattice velocities to physical units [m/s]
        ux_[c]  = ux * lat_to_phys_;
        uy_[c]  = uy * lat_to_phys_;
        uz_[c]  = uz * lat_to_phys_;

        // Store pressure (cell-centered)
        // p = rho * cs2 in lattice units; scale to physical [Pa]
        // For incompressible flow, use dimensionless pressure as relative value
        grid_.p[grid_.c_idx(i, j, k)] = rho * CS2;
    }

    // Convert cell-centered velocities to face-centered (MAC grid)
    // u on x-faces (Nx+1, Ny, Nz)
    for (int i = 0; i <= Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++) {
        double u_face;
        if (i == 0) {
            u_face = ux_[c_idx(0, j, k)];
        } else if (i == Nx) {
            u_face = ux_[c_idx(Nx-1, j, k)];
        } else {
            u_face = 0.5 * (ux_[c_idx(i-1, j, k)] + ux_[c_idx(i, j, k)]);
        }
        grid_.u[grid_.u_idx(i, j, k)] = u_face;
    }

    // v on y-faces (Nx, Ny+1, Nz)
    for (int i = 0; i < Nx; i++)
    for (int j = 0; j <= Ny; j++)
    for (int k = 0; k < Nz; k++) {
        double v_face;
        if (j == 0) {
            v_face = uy_[c_idx(i, 0, k)];
        } else if (j == Ny) {
            v_face = uy_[c_idx(i, Ny-1, k)];
        } else {
            v_face = 0.5 * (uy_[c_idx(i, j-1, k)] + uy_[c_idx(i, j, k)]);
        }
        grid_.v[grid_.v_idx(i, j, k)] = v_face;
    }

    // w on z-faces (Nx, Ny, Nz+1)
    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k <= Nz; k++) {
        double w_face;
        if (k == 0) {
            w_face = uz_[c_idx(i, j, 0)];
        } else if (k == Nz) {
            w_face = uz_[c_idx(i, j, Nz-1)];
        } else {
            w_face = 0.5 * (uz_[c_idx(i, j, k-1)] + uz_[c_idx(i, j, k)]);
        }
        grid_.w[grid_.w_idx(i, j, k)] = w_face;
    }
}

// ---------------------------------------------------------------------------
// compute_divergence_max: uses face-centered velocities written to grid
// ---------------------------------------------------------------------------

double LBMSolver::compute_divergence_max() const {
    double max_div = 0.0;
    for (int i = 0; i < grid_.Nx; i++)
    for (int j = 0; j < grid_.Ny; j++)
    for (int k = 0; k < grid_.Nz; k++) {
        CellType ct = grid_.cell_type[grid_.c_idx(i, j, k)];
        if (ct != CellType::FLUID) continue;
        double div = std::abs(
            (grid_.u[grid_.u_idx(i+1,j,k)] - grid_.u[grid_.u_idx(i,j,k)]) / grid_.dxv[i]
          + (grid_.v[grid_.v_idx(i,j+1,k)] - grid_.v[grid_.v_idx(i,j,k)]) / grid_.dyv[j]
          + (grid_.w[grid_.w_idx(i,j,k+1)] - grid_.w[grid_.w_idx(i,j,k)]) / grid_.dzv[k]);
        max_div = std::max(max_div, div);
    }
    return max_div;
}

// ---------------------------------------------------------------------------
// step: one full LBM iteration
// ---------------------------------------------------------------------------

StepResult LBMSolver::step() {
    // Save previous face-centered velocities for convergence check
    u_prev_ = grid_.u;
    v_prev_ = grid_.v;
    w_prev_ = grid_.w;

    // 1. Set inflow cells to equilibrium (force prescribed velocity)
    apply_inflow_bc();

    // 2. BGK collision
    collision();

    // 3. Streaming + bounce-back
    streaming();

    // 4. Zero-gradient outflow BC (applied after stream)
    apply_outflow_bc();

    // 5. Extract macroscopic fields and write to grid
    extract_macroscopic();

    // Compute velocity change for convergence (velocities are physical [m/s])
    double vel_change = 0.0;
    for (size_t idx = 0; idx < grid_.u.size(); idx++)
        vel_change = std::max(vel_change, std::abs(grid_.u[idx] - u_prev_[idx]));
    for (size_t idx = 0; idx < grid_.v.size(); idx++)
        vel_change = std::max(vel_change, std::abs(grid_.v[idx] - v_prev_[idx]));
    for (size_t idx = 0; idx < grid_.w.size(); idx++)
        vel_change = std::max(vel_change, std::abs(grid_.w[idx] - w_prev_[idx]));

    // Normalize by dt_ (same convention as FluidSolver)
    vel_change /= dt_;

    step_++;

    bool converged = (vel_change < convergence_tol_);

    return {step_, vel_change, 0.0, compute_divergence_max(), converged, 0.0};
}

// ---------------------------------------------------------------------------
// run: loop until convergence or max_steps
// ---------------------------------------------------------------------------

StepResult LBMSolver::run(std::function<void(const StepResult&)> callback) {
    StepResult result{};
    for (int s = 0; s < max_steps_; s++) {
        result = step();
        if (callback) callback(result);
        if (result.converged) break;
    }
    return result;
}

}  // namespace ffroom
