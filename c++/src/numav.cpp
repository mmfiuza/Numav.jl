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
    const uint64_t* const vei_to_ivpg,
    const uint64_t vei_count,
    const Cmplx* const ivpg_to_density,
    const Cmplx* const ivpg_to_soundspeed,
    const uint64_t ivpg_count,
    // surface materials
    const uint64_t* const isei_to_ni,
    const uint64_t* const isei_to_ispgi,
    const uint64_t isei_count,
    const Cmplx* const ispgi_to_impedance,
    const uint64_t ispgi_count,
    // volume velocity
    const uint64_t* const vpi_to_ni,
    const Cmplx* const vpi_to_volvel,
    const uint64_t vpi_count,
    // surface velocity
    const uint64_t* const vsei_to_ni,
    const uint64_t* const vsei_to_ispgv,
    const uint64_t vsei_count,
    const Cmplx* const ispgv_to_velocity,
    const uint64_t ispgv_count,
    // pressure
    const uint64_t* const pni_to_ni,
    const uint64_t pni_count,
    const uint64_t* const pvi_to_pni_count,
    const Cmplx* const pvi_to_pressure,
    const uint64_t pvi_count,
    // export
    const char* const hdf5_file_path,
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
    s._vei_to_ivpg = fz::SafePtr<const uint64_t>::make_view(
        vei_to_ivpg, vei_count
    );
    s._vei_count = vei_count;
    s._ivpg_to_density = fz::SafePtr<const Cmplx>::make_view(
        ivpg_to_density, ivpg_count*fi_count
    );
    s._ivpg_to_soundspeed = fz::SafePtr<const Cmplx>::make_view(
        ivpg_to_soundspeed, ivpg_count*fi_count
    );
    s._ivpg_count = ivpg_count;

    // surface materials
    s._isei_to_ni =
        fz::SafePtr<const std::array<uint64_t, ENIS_COUNT<S,O>>>::make_view(
            reinterpret_cast<const std::array<uint64_t, ENIS_COUNT<S,O>>*>(
                isei_to_ni
            ),
            isei_count*ENIS_COUNT<S,O>
        );
    s._isei_to_ispgi = fz::SafePtr<const uint64_t>::make_view(
        isei_to_ispgi, isei_count
    );
    s._isei_count = isei_count;
    s._ispgi_to_impedance = fz::SafePtr<const Cmplx>::make_view(
        ispgi_to_impedance, ispgi_count*fi_count
    );
    s._ispgi_count = ispgi_count;

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
    s._vsei_to_ispgv = fz::SafePtr<const uint64_t>::make_view(
        vsei_to_ispgv, vsei_count
    );
    s._vsei_count = vsei_count;
    s._ispgv_to_velocity = fz::SafePtr<const Cmplx>::make_view(
        ispgv_to_velocity, ispgv_count*fi_count
    );
    s._ispgv_count = ispgv_count;

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
    s._hdf5_file_path = hdf5_file_path;

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
    const char* const, \
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
