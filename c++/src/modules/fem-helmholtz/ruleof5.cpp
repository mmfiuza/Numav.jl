// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "modules/fem-helmholtz/fem-helmholtz.hpp"

namespace numav {

template<ElementShape S, ElementOrder O>
SimulationFemHelmholtz<S,O>::SimulationFemHelmholtz() {
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
    for (uint64_t ivr = 0UL; ivr != _ivr_count; ++ivr) {
        _ivr_to_stif_fi_part[ivr].free();
        _ivr_to_mass_fi_part[ivr].free();
        _ivr_to_ptr_in_a[ivr].free();
    }
    _ivr_to_stif_fi_part.free();
    _ivr_to_mass_fi_part.free();
    _ivr_to_ptr_in_a.free();
    for (uint64_t isri = 0UL; isri != _isri_count; ++isri) {
        _isri_to_damp_fi_part[isri].free();
        _isri_to_ptr_in_a[isri].free();
    }
    _isri_to_damp_fi_part.free();
    _isri_to_ptr_in_a.free();
    _vpi_to_ptr_in_b.free();
    for (uint64_t isrv = 0UL; isrv != _isrv_count; ++isrv) {
        _isrv_to_forc_fi_part[isrv].free();
        _isrv_to_ptr_in_b[isrv].free();
    }
    _isrv_to_forc_fi_part.free();
    _isrv_to_ptr_in_b.free();
    for (uint64_t pvi = 0UL; pvi != _pvi_count; ++pvi) {
        _pvi_to_ptr_in_a[pvi].free();
        _pvi_to_ptr_in_b[pvi].free();
    }
    _pvi_to_ptr_in_a.free();
    _pvi_to_ptr_in_b.free();
    // _x.free();
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
