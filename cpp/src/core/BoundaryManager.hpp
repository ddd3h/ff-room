#pragma once
#include "Grid.hpp"
#include <array>
#include <vector>

namespace ffroom {

struct FanBC {
    int i_min, i_max;
    int j_min, j_max;
    int k_min, k_max;
    std::array<double, 3> vel;
};

// Outflow patch: axis-aligned face on domain boundary.
// axis: 0=x, 1=y, 2=z   side: 0=low, 1=high
struct OutflowBC {
    int axis;  // 0=x, 1=y, 2=z
    int side;  // 0=low-face, 1=high-face
    int a_min, a_max;  // range in first transverse axis
    int b_min, b_max;  // range in second transverse axis
};

// Opening (window/door): pressure-outflow BC + bidirectional temperature exchange.
// When air flows in from outside: T = T_outside (Dirichlet).
// When air flows out:             T = extrapolated from interior (Neumann).
struct OpeningBC : OutflowBC {
    double T_outside = 20.0;  // outside air temperature [°C]
};

class BoundaryManager {
public:
    void add_solid_box(Grid& grid,
                       int i0, int i1,
                       int j0, int j1,
                       int k0, int k1);

    void add_fan(Grid& grid, const FanBC& fan);

    // Add an outflow patch. Cells on the patch are marked OUTFLOW.
    void add_outflow(Grid& grid, const OutflowBC& obc);

    // Add a window/door opening. Marks cells OUTFLOW + registers temperature BC.
    void add_opening(Grid& grid, const OpeningBC& obc);

    // Apply no-slip on all SOLID faces and domain walls.
    // Skips faces that belong to OUTFLOW cells (those use extrapolation).
    void apply_noslip(Grid& grid);

    void apply_inflow(Grid& grid);

    // Zero-gradient extrapolation on all registered outflow patches.
    void apply_outflow(Grid& grid);

    // Set temperature at opening cells based on current face-velocity direction.
    // Inflow (outside air enters): T = T_outside.
    // Outflow (room air exits):    T = extrapolated from adjacent interior cell.
    void apply_opening_temperature(Grid& grid);

    const std::vector<FanBC>&     fans()     const { return fans_; }
    const std::vector<OutflowBC>& outflows() const { return outflows_; }
    const std::vector<OpeningBC>& openings() const { return openings_; }

private:
    std::vector<FanBC>     fans_;
    std::vector<OutflowBC> outflows_;
    std::vector<OpeningBC> openings_;

    bool _is_outflow_face_x(const Grid& g, int j, int k, bool high) const;
    bool _is_outflow_face_y(const Grid& g, int i, int k, bool high) const;
    bool _is_outflow_face_z(const Grid& g, int i, int j, bool high) const;
};

}  // namespace ffroom
