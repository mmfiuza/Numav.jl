// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include <cstdint>
#include <complex>

namespace numav {

    // type aliases
    using Float = double;
    using Cmplx = typename std::complex<Float>;

    enum class ElementShape : uint64_t {
        TETRAHEDRON
    };

    enum class ElementOrder : uint64_t {
        LINEAR,
        QUADRATIC
    };

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
    );

} // namespace numav
