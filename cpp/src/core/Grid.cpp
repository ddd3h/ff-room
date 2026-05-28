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

    // Initialize uniform face coordinates and spacings
    xs.resize(nx + 1);
    ys.resize(ny + 1);
    zs.resize(nz + 1);
    for (int i = 0; i <= nx; i++) xs[i] = i * dx;
    for (int j = 0; j <= ny; j++) ys[j] = j * dy;
    for (int k = 0; k <= nz; k++) zs[k] = k * dz;

    dxv.assign(nx, dx);
    dyv.assign(ny, dy);
    dzv.assign(nz, dz);
}

void Grid::set_face_coords(std::vector<double> xs_, std::vector<double> ys_,
                           std::vector<double> zs_) {
    xs = std::move(xs_);
    ys = std::move(ys_);
    zs = std::move(zs_);
    for (int i = 0; i < Nx; i++) dxv[i] = xs[i+1] - xs[i];
    for (int j = 0; j < Ny; j++) dyv[j] = ys[j+1] - ys[j];
    for (int k = 0; k < Nz; k++) dzv[k] = zs[k+1] - zs[k];
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
