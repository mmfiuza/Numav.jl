// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "modules/fem-helmholtz/fem-helmholtz.hpp"
#include "common/utils.hpp"

#include <cmath>

namespace numav {

template<ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_clear_data_not_used_in_freq_iterations()
{
    _ni_connections.free();
}

template<ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_solve_systems()
{
    _clear_data_not_used_in_freq_iterations();

    for (uint64_t fi = 0UL; fi != _fi_count; ++fi)
    {
        // frequency calculations
        const Float freq = _fi_to_freq[fi];
        const Float omega = 2_F * PI * freq;
        const Float omega_squared = omega * omega;
        if (freq == 0_F) {
            _x.fill(Cmplx(0_F, 0_F));
            goto write_results_to_file_and_continue;
        }

        // reset values of A matrix and b vector
        _a_vals.fill(Cmplx(0_F, 0_F));
        #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
            _a_diag.fill(Cmplx(0_F, 0_F));
        #endif
        _b_vals.fill(Cmplx(0_F, 0_F));

        // add point volume velocity to b vector
        for (uint64_t vpi = 0UL; vpi != _vpi_count; ++vpi) {
            const Cmplx volvel = _vpi_to_volvel[fi*_vpi_count + vpi];
            *_vpi_to_ptr_in_b[vpi] += Cmplx(0_F, omega) * volvel;
        }

        // add surface velocity to b vector
        for (uint64_t isrv = 0UL; isrv != _isrv_count; ++isrv)
        {
            const Cmplx velocity = _isrv_to_velocity[fi*_isrv_count + isrv];
            const Cmplx fd_part = Cmplx(0_F, omega) * velocity;
            const uint64_t fipi_count = _isrv_to_ptr_in_b[isrv].size();
            for (uint64_t fipi = 0UL; fipi != fipi_count; ++fipi) {
                *_isrv_to_ptr_in_b[isrv][fipi] +=
                    _isrv_to_forc_fi_part[isrv][fipi] * fd_part;
            }
        }

        // add damping matrix to a
        for (uint64_t isri = 0UL; isri != _isri_count; ++isri)
        {
            const Cmplx impedance = _isri_to_impedance[fi*_isri_count+isri];
            const Cmplx damp_fd_part = Cmplx(0_F, omega) / impedance;
            const uint64_t fipi_count = _isri_to_ptr_in_a[isri].size();
            for (uint64_t fipi = 0UL; fipi != fipi_count; ++fipi) {
                *_isri_to_ptr_in_a[isri][fipi] +=
                    _isri_to_damp_fi_part[isri][fipi] * damp_fd_part;
            }
        }

        // add stiffness and mass matrix to a
        for (uint64_t ivr = 0UL; ivr != _ivr_count; ++ivr)
        {
            const Cmplx density = _ivr_to_density[fi*_ivr_count + ivr];
            const Cmplx soundspeed = _ivr_to_soundspeed[fi*_ivr_count+ivr];
            const Cmplx stif_fd_part = Cmplx(1_F, 0_F) / density;
            const Cmplx mass_fd_part =
                omega_squared / (density * soundspeed * soundspeed);
            
            const uint64_t fipi_count = _ivr_to_ptr_in_a[ivr].size();
            for (uint64_t fipi = 0UL; fipi != fipi_count; ++fipi) {
                *_ivr_to_ptr_in_a[ivr][fipi] +=
                    _ivr_to_stif_fi_part[ivr][fipi] * stif_fd_part +
                    _ivr_to_mass_fi_part[ivr][fipi] * mass_fd_part;
            }
        }

        // add pressure to a and b
        for (uint64_t pvi = 0UL; pvi != _pvi_count; ++pvi)
        {
            const Cmplx pressure = _pvi_to_pressure[fi*_pvi_count + pvi];
            const uint64_t fipi_count = _pvi_to_pni_count[pvi];
            for (uint64_t fipi = 0UL; fipi != fipi_count; ++fipi) {
                *_pvi_to_ptr_in_a[pvi][fipi] += PENALTY_METHOD_CONSTANT;
                *_pvi_to_ptr_in_b[pvi][fipi] +=
                    PENALTY_METHOD_CONSTANT * pressure;
            }
        }

        // solve
        #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
            _solve_using_internal_solver();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
            _solve_using_eigen_solver();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
            _solve_using_onemkl_solver();
        #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
            _solve_using_mumps_solver();
        #else
            static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
        #endif

        write_results_to_file_and_continue:

        // external call
        _call_after_every_iteration();
    }

    // terminate solver
    #if NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        _terminate_onemkl_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        _terminate_mumps_solver();
    #endif
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
