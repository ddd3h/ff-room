#include "FanModel.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ffroom {

void DiskFanBC::validate() const {
    for (double x : center)
        if (!std::isfinite(x)) throw std::invalid_argument("Fan center must be finite");
    double norm2 = 0.0;
    for (double x : direction) {
        if (!std::isfinite(x)) throw std::invalid_argument("Fan direction must be finite");
        norm2 += x*x;
    }
    if (!std::isfinite(norm2) || norm2 < 1e-24)
        throw std::invalid_argument("Fan direction must be nonzero");
    if (!std::isfinite(radius) || radius <= 0.0 ||
        !std::isfinite(velocity) || velocity < 0.0 ||
        !std::isfinite(hub_ratio) || hub_ratio < 0.0 || hub_ratio >= 1.0 ||
        !std::isfinite(swirl_ratio))
        throw std::invalid_argument("Invalid fan radius, velocity, hub ratio, or swirl ratio");
}

std::array<double, 3> DiskFanBC::normal() const {
    validate();
    double length = std::hypot(direction[0], direction[1], direction[2]);
    return {direction[0]/length, direction[1]/length, direction[2]/length};
}

double DiskFanBC::radial_distance_squared(const std::array<double, 3>& point) const {
    auto n = normal();
    double s = 0.0;
    for (int a = 0; a < 3; ++a) s += (point[a]-center[a])*n[a];
    double r2 = 0.0;
    for (int a = 0; a < 3; ++a) {
        double radial = point[a]-center[a]-s*n[a];
        r2 += radial*radial;
    }
    return r2;
}

std::array<double, 3> DiskFanBC::velocity_at(const std::array<double, 3>& point) const {
    auto n = normal();
    double s = 0.0;
    for (int a = 0; a < 3; ++a) s += (point[a]-center[a])*n[a];
    std::array<double, 3> r;
    double r2 = 0.0;
    for (int a = 0; a < 3; ++a) {
        r[a] = point[a]-center[a]-s*n[a];
        r2 += r[a]*r[a];
    }
    double outer2 = radius*radius, hub2 = outer2*hub_ratio*hub_ratio;
    if (r2 >= outer2 || r2 < hub2) return {0.0, 0.0, 0.0};
    // q is uniform in annular area; integral_0^1 6q(1-q)dq = 1.
    double q = std::clamp((r2-hub2)/(outer2-hub2), 0.0, 1.0);
    double axial = velocity * (uniform_profile ? 1.0 : 6.0*q*(1.0-q));
    std::array<double, 3> tangent{
        n[1]*r[2]-n[2]*r[1], n[2]*r[0]-n[0]*r[2], n[0]*r[1]-n[1]*r[0]};
    return {axial*(n[0]+swirl_ratio*tangent[0]/radius),
            axial*(n[1]+swirl_ratio*tangent[1]/radius),
            axial*(n[2]+swirl_ratio*tangent[2]/radius)};
}

}  // namespace ffroom
