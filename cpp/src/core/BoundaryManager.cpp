#include "BoundaryManager.hpp"
#include <algorithm>

namespace ffroom {

void BoundaryManager::add_solid_box(Grid& grid,
                                     int i0, int i1,
                                     int j0, int j1,
                                     int k0, int k1) {
    for (int i = std::max(0,i0); i <= std::min(grid.Nx-1,i1); i++)
    for (int j = std::max(0,j0); j <= std::min(grid.Ny-1,j1); j++)
    for (int k = std::max(0,k0); k <= std::min(grid.Nz-1,k1); k++)
        grid.cell_type[grid.c_idx(i,j,k)] = CellType::SOLID;
}

void BoundaryManager::add_fan(Grid& grid, const FanBC& fan) {
    fans_.push_back(fan);
    for (int i = std::max(0,fan.i_min); i <= std::min(grid.Nx-1,fan.i_max); i++)
    for (int j = std::max(0,fan.j_min); j <= std::min(grid.Ny-1,fan.j_max); j++)
    for (int k = std::max(0,fan.k_min); k <= std::min(grid.Nz-1,fan.k_max); k++)
        grid.cell_type[grid.c_idx(i,j,k)] = CellType::INFLOW;
}

void BoundaryManager::add_outflow(Grid& grid, const OutflowBC& obc) {
    outflows_.push_back(obc);
    // Mark cells on the patch as OUTFLOW
    const int Nx = grid.Nx, Ny = grid.Ny, Nz = grid.Nz;
    if (obc.axis == 0) {
        int i = (obc.side == 0) ? 0 : Nx-1;
        for (int j = obc.a_min; j <= obc.a_max; j++)
        for (int k = obc.b_min; k <= obc.b_max; k++)
            grid.cell_type[grid.c_idx(i,j,k)] = CellType::OUTFLOW;
    } else if (obc.axis == 1) {
        int j = (obc.side == 0) ? 0 : Ny-1;
        for (int i = obc.a_min; i <= obc.a_max; i++)
        for (int k = obc.b_min; k <= obc.b_max; k++)
            grid.cell_type[grid.c_idx(i,j,k)] = CellType::OUTFLOW;
    } else {
        int k = (obc.side == 0) ? 0 : Nz-1;
        for (int i = obc.a_min; i <= obc.a_max; i++)
        for (int j = obc.b_min; j <= obc.b_max; j++)
            grid.cell_type[grid.c_idx(i,j,k)] = CellType::OUTFLOW;
    }
}

// Helper: is cell (j,k) on the x-boundary an OUTFLOW cell?
bool BoundaryManager::_is_outflow_face_x(const Grid& g, int j, int k, bool high) const {
    int i = high ? g.Nx-1 : 0;
    return g.cell_type[g.c_idx(i,j,k)] == CellType::OUTFLOW;
}
bool BoundaryManager::_is_outflow_face_y(const Grid& g, int i, int k, bool high) const {
    int j = high ? g.Ny-1 : 0;
    return g.cell_type[g.c_idx(i,j,k)] == CellType::OUTFLOW;
}
bool BoundaryManager::_is_outflow_face_z(const Grid& g, int i, int j, bool high) const {
    int k = high ? g.Nz-1 : 0;
    return g.cell_type[g.c_idx(i,j,k)] == CellType::OUTFLOW;
}

void BoundaryManager::apply_noslip(Grid& grid) {
    // x=0 wall (west)
    for (int j = 0; j < grid.Ny; j++)
    for (int k = 0; k < grid.Nz; k++)
        if (!_is_outflow_face_x(grid, j, k, false))
            grid.u[grid.u_idx(0, j, k)] = 0.0;

    // x=Lx wall (east)
    for (int j = 0; j < grid.Ny; j++)
    for (int k = 0; k < grid.Nz; k++)
        if (!_is_outflow_face_x(grid, j, k, true))
            grid.u[grid.u_idx(grid.Nx, j, k)] = 0.0;

    // y=0 wall (south)
    for (int i = 0; i < grid.Nx; i++)
    for (int k = 0; k < grid.Nz; k++)
        if (!_is_outflow_face_y(grid, i, k, false))
            grid.v[grid.v_idx(i, 0, k)] = 0.0;

    // y=Ly wall (north)
    for (int i = 0; i < grid.Nx; i++)
    for (int k = 0; k < grid.Nz; k++)
        if (!_is_outflow_face_y(grid, i, k, true))
            grid.v[grid.v_idx(i, grid.Ny, k)] = 0.0;

    // z=0 floor
    for (int i = 0; i < grid.Nx; i++)
    for (int j = 0; j < grid.Ny; j++)
        if (!_is_outflow_face_z(grid, i, j, false))
            grid.w[grid.w_idx(i, j, 0)] = 0.0;

    // z=Lz ceiling
    for (int i = 0; i < grid.Nx; i++)
    for (int j = 0; j < grid.Ny; j++)
        if (!_is_outflow_face_z(grid, i, j, true))
            grid.w[grid.w_idx(i, j, grid.Nz)] = 0.0;

    // Internal SOLID faces
    for (int i = 0; i < grid.Nx; i++)
    for (int j = 0; j < grid.Ny; j++)
    for (int k = 0; k < grid.Nz; k++) {
        if (grid.cell_type[grid.c_idx(i,j,k)] != CellType::SOLID) continue;
        grid.u[grid.u_idx(i,   j, k)] = 0.0;
        grid.u[grid.u_idx(i+1, j, k)] = 0.0;
        grid.v[grid.v_idx(i, j,   k)] = 0.0;
        grid.v[grid.v_idx(i, j+1, k)] = 0.0;
        grid.w[grid.w_idx(i, j, k  )] = 0.0;
        grid.w[grid.w_idx(i, j, k+1)] = 0.0;
    }
}

void BoundaryManager::apply_inflow(Grid& grid) {
    for (const auto& fan : fans_) {
        for (int i = std::max(0,fan.i_min); i <= std::min(grid.Nx-1,fan.i_max); i++)
        for (int j = std::max(0,fan.j_min); j <= std::min(grid.Ny-1,fan.j_max); j++)
        for (int k = std::max(0,fan.k_min); k <= std::min(grid.Nz-1,fan.k_max); k++) {
            if (fan.vel[0] != 0.0) {
                grid.u[grid.u_idx(i,   j, k)] = fan.vel[0];
                grid.u[grid.u_idx(i+1, j, k)] = fan.vel[0];
            }
            if (fan.vel[1] != 0.0) {
                grid.v[grid.v_idx(i, j,   k)] = fan.vel[1];
                grid.v[grid.v_idx(i, j+1, k)] = fan.vel[1];
            }
            if (fan.vel[2] != 0.0) {
                grid.w[grid.w_idx(i, j, k  )] = fan.vel[2];
                grid.w[grid.w_idx(i, j, k+1)] = fan.vel[2];
            }
        }
    }
}

void BoundaryManager::add_opening(Grid& grid, const OpeningBC& obc) {
    openings_.push_back(obc);
    add_outflow(grid, static_cast<const OutflowBC&>(obc));
}

void BoundaryManager::apply_opening_temperature(Grid& grid) {
    for (const auto& obc : openings_) {
        const int ax = obc.axis, sd = obc.side;
        const double T_out = obc.T_outside;

        if (ax == 0) {
            int i_wall  = sd ? grid.Nx-1 : 0;
            int i_inner = sd ? grid.Nx-2 : 1;
            int i_face  = sd ? grid.Nx   : 0;
            for (int j = obc.a_min; j <= obc.a_max; j++)
            for (int k = obc.b_min; k <= obc.b_max; k++) {
                double u_f = grid.u[grid.u_idx(i_face, j, k)];
                // sd==0: inflow when u>=0 (outside pushes +x); sd==1: inflow when u<=0
                bool in = sd ? (u_f <= 0.0) : (u_f >= 0.0);
                int c = grid.c_idx(i_wall, j, k);
                grid.T[c] = in ? T_out : grid.T[grid.c_idx(i_inner, j, k)];
            }
        } else if (ax == 1) {
            int j_wall  = sd ? grid.Ny-1 : 0;
            int j_inner = sd ? grid.Ny-2 : 1;
            int j_face  = sd ? grid.Ny   : 0;
            for (int i = obc.a_min; i <= obc.a_max; i++)
            for (int k = obc.b_min; k <= obc.b_max; k++) {
                double v_f = grid.v[grid.v_idx(i, j_face, k)];
                bool in = sd ? (v_f <= 0.0) : (v_f >= 0.0);
                int c = grid.c_idx(i, j_wall, k);
                grid.T[c] = in ? T_out : grid.T[grid.c_idx(i, j_inner, k)];
            }
        } else {
            int k_wall  = sd ? grid.Nz-1 : 0;
            int k_inner = sd ? grid.Nz-2 : 1;
            int k_face  = sd ? grid.Nz   : 0;
            for (int i = obc.a_min; i <= obc.a_max; i++)
            for (int j = obc.b_min; j <= obc.b_max; j++) {
                double w_f = grid.w[grid.w_idx(i, j, k_face)];
                bool in = sd ? (w_f <= 0.0) : (w_f >= 0.0);
                int c = grid.c_idx(i, j, k_wall);
                grid.T[c] = in ? T_out : grid.T[grid.c_idx(i, j, k_inner)];
            }
        }
    }
}

void BoundaryManager::apply_outflow(Grid& grid) {
    for (const auto& obc : outflows_) {
        if (obc.axis == 0) {
            if (obc.side == 1) {  // east face: u[Nx,j,k] = u[Nx-1,j,k]
                for (int j = obc.a_min; j <= obc.a_max; j++)
                for (int k = obc.b_min; k <= obc.b_max; k++)
                    grid.u[grid.u_idx(grid.Nx, j, k)] =
                        grid.u[grid.u_idx(grid.Nx-1, j, k)];
            } else {              // west face: u[0,j,k] = u[1,j,k]
                for (int j = obc.a_min; j <= obc.a_max; j++)
                for (int k = obc.b_min; k <= obc.b_max; k++)
                    grid.u[grid.u_idx(0, j, k)] =
                        grid.u[grid.u_idx(1, j, k)];
            }
        } else if (obc.axis == 1) {
            if (obc.side == 1) {
                for (int i = obc.a_min; i <= obc.a_max; i++)
                for (int k = obc.b_min; k <= obc.b_max; k++)
                    grid.v[grid.v_idx(i, grid.Ny, k)] =
                        grid.v[grid.v_idx(i, grid.Ny-1, k)];
            } else {
                for (int i = obc.a_min; i <= obc.a_max; i++)
                for (int k = obc.b_min; k <= obc.b_max; k++)
                    grid.v[grid.v_idx(i, 0, k)] =
                        grid.v[grid.v_idx(i, 1, k)];
            }
        } else {
            if (obc.side == 1) {
                for (int i = obc.a_min; i <= obc.a_max; i++)
                for (int j = obc.b_min; j <= obc.b_max; j++)
                    grid.w[grid.w_idx(i, j, grid.Nz)] =
                        grid.w[grid.w_idx(i, j, grid.Nz-1)];
            } else {
                for (int i = obc.a_min; i <= obc.a_max; i++)
                for (int j = obc.b_min; j <= obc.b_max; j++)
                    grid.w[grid.w_idx(i, j, 0)] =
                        grid.w[grid.w_idx(i, j, 1)];
            }
        }
    }
}

}  // namespace ffroom
