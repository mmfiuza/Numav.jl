// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "modules/fem-helmholtz/fem-helmholtz.hpp"

namespace numav {

template<ElementShape S, ElementOrder O>
void simulate_fem_helmholtz(
    // freq vector
    const Float* const fi_to_freq,
    const uint64_t fi_count,
    // mesh nodes
    const Float* const ni_to_xyz,
    const uint64_t ni_count,
    // volume materials
    const uint64_t* const vei_to_ni,
    const uint64_t* const vei_to_ivr,
    const uint64_t vei_count,
    const Cmplx* const ivr_to_density,
    const Cmplx* const ivr_to_soundspeed,
    const uint64_t ivr_count,
    // surface materials
    const uint64_t* const isei_to_ni,
    const uint64_t* const isei_to_isri,
    const uint64_t isei_count,
    const Cmplx* const isri_to_impedance,
    const uint64_t isri_count,
    // volume velocity
    const uint64_t* const vpi_to_ni,
    const Cmplx* const vpi_to_volvel,
    const uint64_t vpi_count,
    // surface velocity
    const uint64_t* const vsei_to_ni,
    const uint64_t* const vsei_to_isrv,
    const uint64_t vsei_count,
    const Cmplx* const isrv_to_velocity,
    const uint64_t isrv_count,
    // pressure
    const uint64_t* const pni_to_ni,
    const uint64_t pni_count,
    const uint64_t* const pvi_to_pni_count,
    const Cmplx* const pvi_to_pressure,
    const uint64_t pvi_count,
    // export
    Cmplx* const ni_to_solution,
    // other
    void (*call_after_every_iteration)()
) {
    auto s = SimulationFemHelmholtz<S,O>();

    // freq vector
    s._fi_to_freq = fz::SafePtr<const Float>::make_view(
        fi_to_freq, fi_count
    );
    s._fi_count = fi_count;
    
    // mesh nodes
    s._ni_to_xyz = fz::SafePtr<const std::array<Float, DIM>>::make_view(
        reinterpret_cast<const std::array<Float,DIM>*>(ni_to_xyz), DIM*ni_count
    );
    s._ni_count = ni_count;

    // volume materials
    s._vei_to_ni =
        fz::SafePtr<const std::array<uint64_t, ENIV_COUNT<S,O>>>::make_view(
            reinterpret_cast<const std::array<uint64_t, ENIV_COUNT<S,O>>*>(
                vei_to_ni
            ),
            vei_count*ENIV_COUNT<S,O>
        );
    s._vei_to_ivr = fz::SafePtr<const uint64_t>::make_view(
        vei_to_ivr, vei_count
    );
    s._vei_count = vei_count;
    s._ivr_to_density = fz::SafePtr<const Cmplx>::make_view(
        ivr_to_density, ivr_count*fi_count
    );
    s._ivr_to_soundspeed = fz::SafePtr<const Cmplx>::make_view(
        ivr_to_soundspeed, ivr_count*fi_count
    );
    s._ivr_count = ivr_count;

    // surface materials
    s._isei_to_ni =
        fz::SafePtr<const std::array<uint64_t, ENIS_COUNT<S,O>>>::make_view(
            reinterpret_cast<const std::array<uint64_t, ENIS_COUNT<S,O>>*>(
                isei_to_ni
            ),
            isei_count*ENIS_COUNT<S,O>
        );
    s._isei_to_isri = fz::SafePtr<const uint64_t>::make_view(
        isei_to_isri, isei_count
    );
    s._isei_count = isei_count;
    s._isri_to_impedance = fz::SafePtr<const Cmplx>::make_view(
        isri_to_impedance, isri_count*fi_count
    );
    s._isri_count = isri_count;

    // volume velocity
    s._vpi_to_ni = fz::SafePtr<const uint64_t>::make_view(
        vpi_to_ni, vpi_count
    );
    s._vpi_to_volvel = fz::SafePtr<const Cmplx>::make_view(
        vpi_to_volvel, vpi_count*fi_count
    );
    s._vpi_count = vpi_count;

    // surface velocity
    s._vsei_to_ni =
        fz::SafePtr<const std::array<uint64_t, ENIS_COUNT<S,O>>>::make_view(
            reinterpret_cast<const std::array<uint64_t, ENIS_COUNT<S,O>>*>(
                vsei_to_ni
            ),
            vsei_count*ENIS_COUNT<S,O>
        );
    s._vsei_to_isrv = fz::SafePtr<const uint64_t>::make_view(
        vsei_to_isrv, vsei_count
    );
    s._vsei_count = vsei_count;
    s._isrv_to_velocity = fz::SafePtr<const Cmplx>::make_view(
        isrv_to_velocity, isrv_count*fi_count
    );
    s._isrv_count = isrv_count;

    // pressure
    s._pni_to_ni = fz::SafePtr<const uint64_t>::make_view(
        pni_to_ni, pni_count
    );
    s._pni_count = pni_count;
    s._pvi_to_pni_count = fz::SafePtr<const uint64_t>::make_view(
        pvi_to_pni_count, pvi_count
    );
    s._pvi_to_pressure = fz::SafePtr<const Cmplx>::make_view(
        pvi_to_pressure, pvi_count*fi_count
    );
    s._pvi_count = pvi_count;

    // export
    s._x = fz::SafePtr<Cmplx>::make_view(ni_to_solution, ni_count);

    s._call_after_every_iteration = call_after_every_iteration;

    // run
    s._assemble_freq_independent_parts();
    s._solve_systems();
}

#define SIMULATE_FEM_HELMHOLTZ_PARAMS \
    const Float* const, \
    const uint64_t, \
    const Float* const, \
    const uint64_t, \
    const uint64_t* const, \
    const uint64_t* const, \
    const uint64_t, \
    const Cmplx* const, \
    const Cmplx* const, \
    const uint64_t, \
    const uint64_t* const, \
    const uint64_t* const, \
    const uint64_t, \
    const Cmplx* const, \
    const uint64_t, \
    const uint64_t* const, \
    const Cmplx* const, \
    const uint64_t, \
    const uint64_t* const, \
    const uint64_t* const, \
    const uint64_t, \
    const Cmplx* const, \
    const uint64_t, \
    const uint64_t* const, \
    const uint64_t, \
    const uint64_t* const, \
    const Cmplx* const, \
    const uint64_t, \
    Cmplx* const, \
    void (*)()

// Explicit template instantiations
template
void simulate_fem_helmholtz<ElementShape::TETRAHEDRON,ElementOrder::LINEAR>(
    SIMULATE_FEM_HELMHOLTZ_PARAMS
);
template
void simulate_fem_helmholtz<ElementShape::TETRAHEDRON,ElementOrder::QUADRATIC>(
    SIMULATE_FEM_HELMHOLTZ_PARAMS
);

} // namespace numav
