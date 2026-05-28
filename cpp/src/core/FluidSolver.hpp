#pragma once
#include "Grid.hpp"
#include "PoissonSolver.hpp"
#include "BoundaryManager.hpp"
#include <functional>

namespace ffroom {

enum class AdvectionScheme { UPWIND1 = 0, QUICK = 1, LAX_WENDROFF = 2 };

struct FluidSolverParams {
    double rho     = 1.2;     // air density [kg/m3]
    double nu      = 1.5e-5;  // kinematic viscosity [m2/s] (air at 20°C)
    double dt      = 0.01;    // time step [s]
    int    max_steps    = 500;
    double convergence_tol = 1.0;  // max |u_new - u_old| / dt for steady-state

    PoissonSolverParams poisson;

    // Advection scheme (default UPWIND1 preserves existing behavior)
    AdvectionScheme advection = AdvectionScheme::UPWIND1;

    // OpenMP parallelization (default false preserves existing behavior)
    bool use_openmp = false;

    // Use multigrid V-cycle solver instead of CG for pressure Poisson
    bool use_multigrid = false;

    // Thermal simulation
    bool   thermal     = false;       // enable temperature advection-diffusion
    double T_initial   = 20.0;        // initial / reference temperature [°C]
    double T_target    = 1e30;        // stop when mean(T) reaches this; 1e30 = never
    bool   buoyancy    = false;       // Boussinesq buoyancy (hot air rises in +z)
    double g_accel     = 9.81;        // gravitational acceleration [m/s²]
    double beta        = 1.0/293.0;   // thermal expansion coefficient [1/K]
    double k_thermal   = 0.026;       // thermal conductivity of air [W/m/K]
    double cp          = 1005.0;      // specific heat of air [J/kg/K]
};

struct StepResult {
    int    step;
    double velocity_change;   // max |u^{n+1} - u^n| / dt
    double pressure_residual; // Poisson solver residual
    double divergence_max;    // max |∇·u| after correction
    bool   converged;
    double T_mean = 0.0;      // mean temperature of FLUID cells [°C]
};

class FluidSolver {
public:
    FluidSolver(Grid& grid, BoundaryManager& bm, FluidSolverParams params = {});

    // Run one projection-method time step. Returns diagnostics.
    StepResult step();

    // Run until convergence or max_steps. Calls callback each step (may be null).
    StepResult run(std::function<void(const StepResult&)> callback = nullptr);

    double compute_divergence_max() const;
    double compute_T_mean() const;

    const FluidSolverParams& params() const { return params_; }
    int current_step() const { return step_; }

private:
    Grid&            grid_;
    BoundaryManager& bm_;
    FluidSolverParams params_;
    PoissonSolver             poisson_;
    MultigridPoissonSolver    multigrid_poisson_;
    int               step_ = 0;

    std::vector<double> u_prev_, v_prev_, w_prev_;

    void compute_intermediate();
    void build_poisson_rhs(std::vector<double>& rhs) const;
    void apply_correction();

    // Thermal methods
    void apply_buoyancy();      // add g*β*(T-T_ref)*dt to w* (z-velocity)
    void advect_temperature();  // explicit upwind advection + diffusion of T
};

}  // namespace ffroom
