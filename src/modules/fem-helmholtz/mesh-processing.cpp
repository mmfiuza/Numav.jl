// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "common/exception.hpp"
#include "common/log.hpp"
#include "common/hash-functions.hpp"
#include "common/maths.hpp"
#include "common/utils.hpp"

#include <tuple>
#include <fstream>
#include <limits>

namespace numav {

template <ElementOrder O>
uint64_t SimulationFemHelmTet<O>::_get_closest_point(
    const std::array<Float,DIM> point_coords
) {
    Float minimum_distance_squared = std::numeric_limits<Float>::max();
    uint64_t ni_closest = std::numeric_limits<uint64_t>::max();
    for (uint64_t ni = 0UL; ni != _ni_count; ++ni) {
        Float distance_squared = 
            power<2UL>(_ni_to_xyz[ni][0UL] - point_coords[0UL]) +
            power<2UL>(_ni_to_xyz[ni][1UL] - point_coords[1UL]) +
            power<2UL>(_ni_to_xyz[ni][2UL] - point_coords[2UL])
        ;
        if (distance_squared < minimum_distance_squared) {
            minimum_distance_squared = distance_squared;
            ni_closest = ni;
        }
    }
    return ni_closest;
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
