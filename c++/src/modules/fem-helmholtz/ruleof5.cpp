// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "modules/fem-helmholtz/fem-helmholtz.hpp"
#include "common/log.hpp"

namespace numav {

template<ElementShape S, ElementOrder O>
SimulationFemHelmholtz<S,O>::SimulationFemHelmholtz() {
    log::set_level();
    log::set_pattern();
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        _solver = std::make_unique<Eigen::SparseLU<
            Eigen::SparseMatrix<Cmplx, Eigen::ColMajor, Eigen::Index>,
            Eigen::COLAMDOrdering<Eigen::Index>
        >>();
    #endif
}

template<ElementShape S, ElementOrder O>
SimulationFemHelmholtz<S,O>::~SimulationFemHelmholtz() {
    _a_vals.free();
    _b_row_idx.free();
    _b_vals.free();
    for (uint64_t ivpg = 0UL; ivpg != _ivpg_count; ++ivpg) {
        _ivpg_to_stif_fi_part[ivpg].free();
        _ivpg_to_mass_fi_part[ivpg].free();
        _ivpg_to_ptr_in_a[ivpg].free();
    }
    _ivpg_to_stif_fi_part.free();
    _ivpg_to_mass_fi_part.free();
    _ivpg_to_ptr_in_a.free();
    for (uint64_t ispgi = 0UL; ispgi != _ispgi_count; ++ispgi) {
        _ispgi_to_damp_fi_part[ispgi].free();
        _ispgi_to_ptr_in_a[ispgi].free();
    }
    _ispgi_to_damp_fi_part.free();
    _ispgi_to_ptr_in_a.free();
    _vpi_to_ptr_in_b.free();
    for (uint64_t ispgv = 0UL; ispgv != _ispgv_count; ++ispgv) {
        _ispgv_to_forc_fi_part[ispgv].free();
        _ispgv_to_ptr_in_b[ispgv].free();
    }
    _ispgv_to_forc_fi_part.free();
    _ispgv_to_ptr_in_b.free();
    for (uint64_t pvi = 0UL; pvi != _pvi_count; ++pvi) {
        _pvi_to_ptr_in_a[pvi].free();
        _pvi_to_ptr_in_b[pvi].free();
    }
    _pvi_to_ptr_in_a.free();
    _pvi_to_ptr_in_b.free();
    _x.free();
    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        _b_dense.free();
        _a_diag.free();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        _a_row_idx.free();
        _a_col_idx.free();
        _b_row_idx_signed.free();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        _b_dense.free();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        _b_dense.free();
    #endif
}

template<ElementShape S, ElementOrder O>
SimulationFemHelmholtz<S,O>::SimulationFemHelmholtz(
    SimulationFemHelmholtz&& other
) noexcept = default;

template<ElementShape S, ElementOrder O>
SimulationFemHelmholtz<S,O>& SimulationFemHelmholtz<S,O>::operator=(
    SimulationFemHelmholtz&& other
) noexcept = default;

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
