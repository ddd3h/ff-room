#pragma once
#include <array>

namespace ffroom {

// A prescribed outlet profile on a circular plane, not a resolved blade model.
struct DiskFanBC {
    std::array<double, 3> center{0.5, 0.5, 0.5};
    std::array<double, 3> direction{1.0, 0.0, 0.0};
    double velocity = 3.0;       // area-mean axial speed over the blade annulus
    double radius = 0.15;
    double hub_ratio = 0.2;      // hub radius / outer radius
    double swirl_ratio = 0.0;    // signed tangential / axial speed at the rim
    bool uniform_profile = false;

    void validate() const;
    std::array<double, 3> normal() const;
    double radial_distance_squared(const std::array<double, 3>& point) const;
    std::array<double, 3> velocity_at(const std::array<double, 3>& point) const;
};

struct FanInflowCell {
    std::array<int, 3> cell;
    std::array<double, 3> velocity;
};

}  // namespace ffroom
