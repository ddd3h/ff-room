#pragma once
#include "Grid.hpp"
#include "BoundaryManager.hpp"
#include "FluidSolver.hpp"   // for StepResult
#include <functional>
#include <vector>

namespace ffroom {

// D3Q19 Lattice Boltzmann Method solver.
// Uses BGK collision with a single relaxation time tau.
// Isothermal incompressible flow (cs2 = 1/3 in lattice units).
//
// Physical-to-lattice unit conversion:
//   u_lat = u_phys * dt / dx   (lattice velocity, must satisfy Ma << 1)
//   When extracting macroscopic fields, we convert back:
//   u_phys = u_lat * dx / dt
//
// Distribution layout: f[q * (Nx*Ny*Nz) + i*(Ny*Nz) + j*Nz + k]
//   q = 0..18, i = 0..Nx-1, j = 0..Ny-1, k = 0..Nz-1
class LBMSolver {
public:
    // tau:              BGK relaxation time (0.5 < tau < 2.0)
    // dt:               physical time step [s] — used for physical/lattice conversion
    // max_steps, convergence_tol: same semantics as FluidSolver
    LBMSolver(Grid& grid, BoundaryManager& bm, double tau,
              int max_steps, double convergence_tol, double dt = 0.01);

    // Run one LBM time step.  Returns diagnostics.
    StepResult step();

    // Run until convergence or max_steps. Calls callback each step (may be null).
    StepResult run(std::function<void(const StepResult&)> callback = nullptr);

    double compute_divergence_max() const;
    double compute_T_mean() const { return 0.0; }  // isothermal — no thermal field

    int current_step() const { return step_; }

private:
    Grid&            grid_;
    BoundaryManager& bm_;
    double           tau_;
    int              max_steps_;
    double           convergence_tol_;
    double           dt_;       // physical time step [s]
    int              step_ = 0;

    // Velocity scaling: lat = phys * (dt / dx_ref)
    // We use grid_.dx as reference length scale.
    // lat_to_phys = dx_ref / dt; phys_to_lat = dt / dx_ref
    double lat_to_phys_;   // multiply lattice velocity by this to get physical [m/s]
    double phys_to_lat_;   // multiply physical velocity by this to get lattice units

    // Distribution functions (double-buffered)
    std::vector<double> f_;     // current distributions
    std::vector<double> f_tmp_; // post-stream buffer

    // Macroscopic cell-centered velocities (updated each step)
    std::vector<double> ux_, uy_, uz_; // size Nx*Ny*Nz each
    std::vector<double> rho_;           // size Nx*Ny*Nz

    // Previous step velocities for convergence check (stored as face-centered
    // copies of grid_.u/v/w before this step starts)
    std::vector<double> u_prev_, v_prev_, w_prev_;

    // D3Q19 constants
    static constexpr int Q = 19;
    static constexpr double CS2 = 1.0 / 3.0;  // speed of sound squared

    static const int EX[19];
    static const int EY[19];
    static const int EZ[19];
    static const double W[19];
    static const int OPP[19];

    inline int f_idx(int q, int i, int j, int k) const {
        return q * (grid_.Nx * grid_.Ny * grid_.Nz) +
               i * (grid_.Ny * grid_.Nz) +
               j * grid_.Nz +
               k;
    }

    inline int c_idx(int i, int j, int k) const {
        return i * (grid_.Ny * grid_.Nz) + j * grid_.Nz + k;
    }

    // Compute equilibrium distribution for velocity (ux,uy,uz) and density rho
    inline double feq(int q, double rho, double ux, double uy, double uz) const {
        double eu = EX[q]*ux + EY[q]*uy + EZ[q]*uz;
        double usq = ux*ux + uy*uy + uz*uz;
        return W[q] * rho * (1.0 + eu/CS2 + eu*eu/(2.0*CS2*CS2) - usq/(2.0*CS2));
    }

    void initialize_f();           // set f to equilibrium with zero velocity
    void collision();              // BGK collision in-place on f_
    void streaming();              // stream f_ -> f_tmp_, apply bounce-back BCs
    void apply_inflow_bc();        // set inflow cells to eq each step (before collision)
    void apply_outflow_bc();       // zero-gradient copy for outflow cells (after stream)
    void extract_macroscopic();    // compute rho,ux,uy,uz from f_; write to grid

    // Check if cell (i,j,k) is inside domain and not SOLID/OUTFLOW
    inline bool is_fluid_or_inflow(int i, int j, int k) const {
        if (i < 0 || i >= grid_.Nx || j < 0 || j >= grid_.Ny || k < 0 || k >= grid_.Nz)
            return false;
        CellType ct = grid_.cell_type[grid_.c_idx(i,j,k)];
        return ct == CellType::FLUID || ct == CellType::INFLOW;
    }
};

}  // namespace ffroom
