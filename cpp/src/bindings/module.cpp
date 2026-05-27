#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include "../core/Grid.hpp"
#include "../core/PoissonSolver.hpp"
#include "../core/BoundaryManager.hpp"
#include "../core/FluidSolver.hpp"

namespace py = pybind11;
using namespace ffroom;

PYBIND11_MODULE(_ffroom_core, m) {
    m.doc() = "ff-room CFD core (projection method on MAC grid)";

    // CellType enum
    py::enum_<CellType>(m, "CellType")
        .value("FLUID",   CellType::FLUID)
        .value("SOLID",   CellType::SOLID)
        .value("INFLOW",  CellType::INFLOW)
        .value("OUTFLOW", CellType::OUTFLOW)
        .export_values();

    // Grid
    py::class_<Grid>(m, "Grid")
        .def(py::init<int,int,int,double,double,double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"),
             py::arg("lx"), py::arg("ly"), py::arg("lz"))
        .def_readonly("Nx", &Grid::Nx)
        .def_readonly("Ny", &Grid::Ny)
        .def_readonly("Nz", &Grid::Nz)
        .def_readonly("dx", &Grid::dx)
        .def_readonly("dy", &Grid::dy)
        .def_readonly("dz", &Grid::dz)
        // Expose velocity/pressure as numpy arrays (zero-copy view)
        .def_property_readonly("u", [](Grid& g) {
            return py::array_t<double>(
                {(g.Nx+1)*g.Ny*g.Nz}, g.u.data(), py::cast(g));
        })
        .def_property_readonly("v", [](Grid& g) {
            return py::array_t<double>(
                {g.Nx*(g.Ny+1)*g.Nz}, g.v.data(), py::cast(g));
        })
        .def_property_readonly("w", [](Grid& g) {
            return py::array_t<double>(
                {g.Nx*g.Ny*(g.Nz+1)}, g.w.data(), py::cast(g));
        })
        .def_property_readonly("p", [](Grid& g) {
            return py::array_t<double>(
                {g.Nx*g.Ny*g.Nz}, g.p.data(), py::cast(g));
        })
        // Shaped views for easier Python indexing
        .def("u_shaped", [](Grid& g) {
            return py::array_t<double>(
                {g.Nx+1, g.Ny, g.Nz}, g.u.data(), py::cast(g));
        })
        .def("v_shaped", [](Grid& g) {
            return py::array_t<double>(
                {g.Nx, g.Ny+1, g.Nz}, g.v.data(), py::cast(g));
        })
        .def("w_shaped", [](Grid& g) {
            return py::array_t<double>(
                {g.Nx, g.Ny, g.Nz+1}, g.w.data(), py::cast(g));
        })
        .def("p_shaped", [](Grid& g) {
            return py::array_t<double>(
                {g.Nx, g.Ny, g.Nz}, g.p.data(), py::cast(g));
        })
        .def("T_shaped", [](Grid& g) {
            return py::array_t<double>({g.Nx, g.Ny, g.Nz}, g.T.data(), py::cast(g));
        })
        .def("cell_type_shaped", [](Grid& g) {
            // Return as uint8 array
            std::vector<uint8_t> ct(g.cell_type.size());
            for (size_t i = 0; i < ct.size(); i++)
                ct[i] = static_cast<uint8_t>(g.cell_type[i]);
            auto arr = py::array_t<uint8_t>({g.Nx, g.Ny, g.Nz});
            std::memcpy(arr.mutable_data(), ct.data(), ct.size());
            return arr;
        })
        .def("zero_velocity", &Grid::zero_velocity)
        .def("zero_pressure", &Grid::zero_pressure);

    // FanBC
    py::class_<FanBC>(m, "FanBC")
        .def(py::init<>())
        .def_readwrite("i_min", &FanBC::i_min)
        .def_readwrite("i_max", &FanBC::i_max)
        .def_readwrite("j_min", &FanBC::j_min)
        .def_readwrite("j_max", &FanBC::j_max)
        .def_readwrite("k_min", &FanBC::k_min)
        .def_readwrite("k_max", &FanBC::k_max)
        .def_readwrite("vel",   &FanBC::vel);

    // OutflowBC
    py::class_<OutflowBC>(m, "OutflowBC")
        .def(py::init<>())
        .def_readwrite("axis",  &OutflowBC::axis)
        .def_readwrite("side",  &OutflowBC::side)
        .def_readwrite("a_min", &OutflowBC::a_min)
        .def_readwrite("a_max", &OutflowBC::a_max)
        .def_readwrite("b_min", &OutflowBC::b_min)
        .def_readwrite("b_max", &OutflowBC::b_max);

    // OpeningBC (window / door)
    py::class_<OpeningBC, OutflowBC>(m, "OpeningBC")
        .def(py::init<>())
        .def_readwrite("axis",      &OpeningBC::axis)
        .def_readwrite("side",      &OpeningBC::side)
        .def_readwrite("a_min",     &OpeningBC::a_min)
        .def_readwrite("a_max",     &OpeningBC::a_max)
        .def_readwrite("b_min",     &OpeningBC::b_min)
        .def_readwrite("b_max",     &OpeningBC::b_max)
        .def_readwrite("T_outside", &OpeningBC::T_outside);

    // BoundaryManager
    py::class_<BoundaryManager>(m, "BoundaryManager")
        .def(py::init<>())
        .def("add_solid_box",              &BoundaryManager::add_solid_box)
        .def("add_fan",                    &BoundaryManager::add_fan)
        .def("add_outflow",                &BoundaryManager::add_outflow)
        .def("add_opening",                &BoundaryManager::add_opening)
        .def("apply_noslip",               &BoundaryManager::apply_noslip)
        .def("apply_inflow",               &BoundaryManager::apply_inflow)
        .def("apply_outflow",              &BoundaryManager::apply_outflow)
        .def("apply_opening_temperature",  &BoundaryManager::apply_opening_temperature);

    // PoissonSolverParams
    py::class_<PoissonSolverParams>(m, "PoissonSolverParams")
        .def(py::init<>())
        .def_readwrite("max_iter", &PoissonSolverParams::max_iter)
        .def_readwrite("tol",      &PoissonSolverParams::tol);

    // FluidSolverParams
    py::class_<FluidSolverParams>(m, "FluidSolverParams")
        .def(py::init<>())
        .def_readwrite("rho",              &FluidSolverParams::rho)
        .def_readwrite("nu",               &FluidSolverParams::nu)
        .def_readwrite("dt",               &FluidSolverParams::dt)
        .def_readwrite("max_steps",        &FluidSolverParams::max_steps)
        .def_readwrite("convergence_tol",  &FluidSolverParams::convergence_tol)
        .def_readwrite("poisson",          &FluidSolverParams::poisson)
        // Thermal
        .def_readwrite("thermal",          &FluidSolverParams::thermal)
        .def_readwrite("T_initial",        &FluidSolverParams::T_initial)
        .def_readwrite("T_target",         &FluidSolverParams::T_target)
        .def_readwrite("buoyancy",         &FluidSolverParams::buoyancy)
        .def_readwrite("g_accel",          &FluidSolverParams::g_accel)
        .def_readwrite("beta",             &FluidSolverParams::beta)
        .def_readwrite("k_thermal",        &FluidSolverParams::k_thermal)
        .def_readwrite("cp",               &FluidSolverParams::cp);

    // StepResult
    py::class_<StepResult>(m, "StepResult")
        .def_readonly("step",              &StepResult::step)
        .def_readonly("velocity_change",   &StepResult::velocity_change)
        .def_readonly("pressure_residual", &StepResult::pressure_residual)
        .def_readonly("divergence_max",    &StepResult::divergence_max)
        .def_readonly("converged",         &StepResult::converged)
        .def_readonly("T_mean",            &StepResult::T_mean)
        .def("__repr__", [](const StepResult& r) {
            return "StepResult(step=" + std::to_string(r.step)
                + ", vel_change=" + std::to_string(r.velocity_change)
                + ", p_res=" + std::to_string(r.pressure_residual)
                + ", div_max=" + std::to_string(r.divergence_max)
                + ", T_mean=" + std::to_string(r.T_mean)
                + ", converged=" + (r.converged ? "True" : "False") + ")";
        });

    // FluidSolver
    py::class_<FluidSolver>(m, "FluidSolver")
        .def(py::init<Grid&, BoundaryManager&, FluidSolverParams>(),
             py::arg("grid"), py::arg("bm"),
             py::arg("params") = FluidSolverParams{})
        .def("step",            &FluidSolver::step)
        .def("run",             &FluidSolver::run,
             py::arg("callback") = nullptr)
        .def("compute_divergence_max", &FluidSolver::compute_divergence_max)
        .def_property_readonly("current_step", &FluidSolver::current_step);
}
