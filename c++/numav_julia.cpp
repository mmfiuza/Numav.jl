// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "numav/numav.hpp"

namespace {

template<numav::ElementShape S, numav::ElementOrder O>
void wrapped_simulate_fem_helmholtz(
    // freq vector
    const numav::Float* const fi_to_freq,
    const uint64_t fi_count,
    // mesh nodes
    const numav::Float* const ni_to_xyz,
    const uint64_t ni_count,
    // volume materials
    const uint64_t* const vei_to_ni,
    const uint64_t* const vei_to_ivpg,
    const uint64_t vei_count,
    const numav::Cmplx* const ivpg_to_density,
    const numav::Cmplx* const ivpg_to_soundspeed,
    const uint64_t ivpg_count,
    // surface materials
    const uint64_t* const isei_to_ni,
    const uint64_t* const isei_to_ispgi,
    const uint64_t isei_count,
    const numav::Cmplx* const ispgi_to_impedance,
    const uint64_t ispgi_count,
    // volume velocity
    const uint64_t* const vpi_to_ni,
    const numav::Cmplx* const vpi_to_volvel,
    const uint64_t vpi_count,
    // surface velocity
    const uint64_t* const vsei_to_ni,
    const uint64_t* const vsei_to_ispgv,
    const uint64_t vsei_count,
    const numav::Cmplx* const ispgv_to_velocity,
    const uint64_t ispgv_count,
    // pressure
    const uint64_t* const pni_to_ni,
    const uint64_t pni_count,
    const uint64_t* const pvi_to_pni_count,
    const numav::Cmplx* const pvi_to_pressure,
    const uint64_t pvi_count,
    // export
    numav::Cmplx* const ni_to_solution,
    // other
    jlcxx::SafeCFunction call_after_every_iteration_safe
) {
    void (*call_after_every_iteration)() =
        jlcxx::make_function_pointer<void()>(call_after_every_iteration_safe);

    numav::simulate_fem_helmholtz<S, O>(
        fi_to_freq, fi_count,
        ni_to_xyz, ni_count,
        vei_to_ni, vei_to_ivpg, vei_count,
        ivpg_to_density, ivpg_to_soundspeed, ivpg_count,
        isei_to_ni, isei_to_ispgi, isei_count,
        ispgi_to_impedance, ispgi_count,
        vpi_to_ni, vpi_to_volvel, vpi_count,
        vsei_to_ni, vsei_to_ispgv, vsei_count,
        ispgv_to_velocity, ispgv_count,
        pni_to_ni, pni_count,
        pvi_to_pni_count, pvi_to_pressure, pvi_count,
        ni_to_solution,
        call_after_every_iteration
    );
}

} // anonymous namespace

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.method("_cpp_simulate_fem_helmholtz_tetrahedron_linear",
        &wrapped_simulate_fem_helmholtz<
            numav::ElementShape::TETRAHEDRON,
            numav::ElementOrder::LINEAR
        >
    );
    mod.method("_cpp_simulate_fem_helmholtz_tetrahedron_quadratic",
        &wrapped_simulate_fem_helmholtz<
            numav::ElementShape::TETRAHEDRON,
            numav::ElementOrder::QUADRATIC
        >
    );
}
