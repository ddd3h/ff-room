#include "Grid.hpp"
#include <stdexcept>

namespace ffroom {

Grid::Grid(int nx, int ny, int nz, double lx, double ly, double lz)
    : Nx(nx), Ny(ny), Nz(nz),
      dx(lx / nx), dy(ly / ny), dz(lz / nz)
{
    if (nx <= 0 || ny <= 0 || nz <= 0)
        throw std::invalid_argument("Grid dimensions must be positive");

    u.assign((nx+1)*ny*nz, 0.0);
    v.assign(nx*(ny+1)*nz, 0.0);
    w.assign(nx*ny*(nz+1), 0.0);
    p.assign(nx*ny*nz, 0.0);

    u_tmp.assign((nx+1)*ny*nz, 0.0);
    v_tmp.assign(nx*(ny+1)*nz, 0.0);
    w_tmp.assign(nx*ny*(nz+1), 0.0);

    cell_type.assign(nx*ny*nz, CellType::FLUID);
    T.assign(nx*ny*nz, 0.0);
    T_tmp.assign(nx*ny*nz, 0.0);
}

void Grid::zero_velocity() {
    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(w.begin(), w.end(), 0.0);
}

void Grid::zero_pressure() {
    std::fill(p.begin(), p.end(), 0.0);
}

}  // namespace ffroom
