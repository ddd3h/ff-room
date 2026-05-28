#pragma once
#include <vector>
#include <array>
#include <cstdint>

namespace ffroom {

// Cell type flags
enum class CellType : uint8_t {
    FLUID  = 0,
    SOLID  = 1,   // no-slip wall / obstacle
    INFLOW = 2,   // fan face
    OUTFLOW = 3,  // passive outflow (pressure = 0)
};

// MAC (Marker-And-Cell) staggered grid
//   u: face-centered on x-faces  -> (Nx+1, Ny, Nz)
//   v: face-centered on y-faces  -> (Nx, Ny+1, Nz)
//   w: face-centered on z-faces  -> (Nx, Ny, Nz+1)
//   p: cell-centered             -> (Nx, Ny, Nz)
//   cell_type: cell-centered     -> (Nx, Ny, Nz)
//
// Coordinate: x=width(east), y=depth(north), z=height(up)
// Unit: SI (m, m/s, Pa, kg/m3)

struct Grid {
    int Nx, Ny, Nz;
    double dx, dy, dz;  // mean cell size [m] (= Lx/Nx etc.; kept for compat)

    // Variable cell spacings (non-uniform grid)
    std::vector<double> dxv;  // dxv[i] = x_face[i+1] - x_face[i], size Nx
    std::vector<double> dyv;  // dyv[j] = y_face[j+1] - y_face[j], size Ny
    std::vector<double> dzv;  // dzv[k] = z_face[k+1] - z_face[k], size Nz
    // Face coordinates
    std::vector<double> xs;   // x face positions, size Nx+1
    std::vector<double> ys;   // y face positions, size Ny+1
    std::vector<double> zs;   // z face positions, size Nz+1

    // Initialize/recompute spacing from face coordinate arrays
    void set_face_coords(std::vector<double> xs_, std::vector<double> ys_,
                         std::vector<double> zs_);

    std::vector<double> u;  // x-velocity on x-faces
    std::vector<double> v;  // y-velocity on y-faces
    std::vector<double> w;  // z-velocity on z-faces
    std::vector<double> p;  // pressure, cell-centered

    std::vector<double> u_tmp, v_tmp, w_tmp;  // intermediate velocity u*

    std::vector<CellType> cell_type;

    std::vector<double> T;      // temperature, cell-centered [°C]
    std::vector<double> T_tmp;  // scratch buffer for temperature advection

    Grid(int nx, int ny, int nz, double lx, double ly, double lz);

    // Flat index helpers
    // u: (Nx+1)*Ny*Nz
    inline int u_idx(int i, int j, int k) const { return i*(Ny*Nz) + j*Nz + k; }
    // v: Nx*(Ny+1)*Nz
    inline int v_idx(int i, int j, int k) const { return i*((Ny+1)*Nz) + j*Nz + k; }
    // w: Nx*Ny*(Nz+1)
    inline int w_idx(int i, int j, int k) const { return i*(Ny*(Nz+1)) + j*(Nz+1) + k; }
    // p / cell_type: Nx*Ny*Nz
    inline int c_idx(int i, int j, int k) const { return i*(Ny*Nz) + j*Nz + k; }

    void zero_velocity();
    void zero_pressure();
};

}  // namespace ffroom
