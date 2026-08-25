// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/modules/fem-helmholtz/compile-options.hpp"

#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <memory>
#include <optional>

#include "Eigen/Eigen"
#include "Eigen/OrderingMethods"
#include "SafePtr.hpp"
#include "indicators/progress_bar.hpp"
#include "H5Cpp.h"
#include "boost/bimap.hpp"
#include "common/log.hpp"

namespace numav {

constexpr uint64_t DIM = 3UL; // dimension count in space
constexpr Cmplx PENALTY_METHOD_CONSTANT = Cmplx(1e12_F, 0_F);

enum class Equation : uint64_t {
    HELMHOLTZ
};

enum class ElementShape : uint64_t {
    TETRAHEDRON
};

enum class ElementOrder : uint64_t {
    LINEAR,
    QUADRATIC
};

template<ElementOrder O> constexpr uint64_t ENIS_COUNT = [] {
    if constexpr (O == ElementOrder::LINEAR) { return 3UL; }
    if constexpr (O == ElementOrder::QUADRATIC) { return 6UL; }
    return 0UL;
}();

template<ElementOrder O> constexpr uint64_t EXTRA_ENIS_COUNT = [] {
    return ENIS_COUNT<O> - ENIS_COUNT<ElementOrder::LINEAR>;
}();

template<ElementOrder O> constexpr uint64_t ENIV_COUNT = [] {
    if constexpr (O == ElementOrder::LINEAR) { return 4UL;  }
    if constexpr (O == ElementOrder::QUADRATIC) { return 10UL; }
    return 0UL;
}();

template<ElementOrder O> constexpr uint64_t EXTRA_ENIV_COUNT = [] {
    return ENIV_COUNT<O> - ENIV_COUNT<ElementOrder::LINEAR>;
}();

template<ElementOrder O>
class Simulation<
    NumericalMethod::FEM,
    Equation::HELMHOLTZ,
    ElementShape::TETRAHEDRON,
    O
> {
public:
    Simulation();
    ~Simulation();
    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept;
    Simulation& operator=(Simulation&&) noexcept;

    void _allocate_a();
    void _allocate_b();
    void _allocate_x();
    void _assemble_fi_part_for_point_velocity();
    void _assemble_fi_part_for_sfc_velocity();
    void _assemble_fi_part_for_sfc_impedance();
    void _assemble_fi_part_for_vol_elements();
    void _assemble_fi_part_for_pressure();
    void _assemble_freq_independent_parts();
    void _solve_systems();
    H5::DataSet _begin_hdf5_file();
    void _write_solution_for_one_freq(H5::DataSet& ds, const uint64_t fi);
    void _clear_data_not_used_in_freq_iterations();
    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        void _define_sparsity_pattern_using_internal_solver();
        void _solve_using_internal_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        void _define_sparsity_pattern_using_eigen_solver();
        void _solve_using_eigen_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        void _define_sparsity_pattern_using_onemkl_solver();
        void _solve_using_onemkl_solver();
        void _terminate_onemkl_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        void _define_sparsity_pattern_using_mumps_solver();
        void _solve_using_mumps_solver();
        void _terminate_mumps_solver();
    #endif
    
    H5::H5File _hdf5_file;

    std::string _hdf5_file_path;
    
    fz::SafePtr<const Float> _fi_to_freq;

    fz::SafePtr<const std::array<Float, DIM>> _ni_to_xyz;

    fz::SafePtr<const std::array<uint64_t, ENIV_COUNT<O>>> _vei_to_ni;
    fz::SafePtr<const uint64_t> _vei_to_ivpg;
    fz::SafePtr<const Cmplx> _ivpg_to_density;
    fz::SafePtr<const Cmplx> _ivpg_to_soundspeed;

    fz::SafePtr<const std::array<uint64_t, ENIS_COUNT<O>>> _isei_to_ni;
    fz::SafePtr<const uint64_t> _isei_to_ispgi;
    fz::SafePtr<const Cmplx> _ispgi_to_impedance;

    fz::SafePtr<const uint64_t> _vpi_to_ni;
    fz::SafePtr<const Cmplx> _vpi_to_volvel;

    fz::SafePtr<const std::array<uint64_t, ENIS_COUNT<O>>> _vsei_to_ni;
    fz::SafePtr<const uint64_t> _vsei_to_ispgv;
    fz::SafePtr<const Cmplx> _ispgv_to_velocity;
    
    fz::SafePtr<const uint64_t> _pni_to_ni;
    fz::SafePtr<const uint64_t> _pvi_to_pni_count;
    fz::SafePtr<const Cmplx> _pvi_to_pressure;

    fz::SafePtr<std::pair<uint64_t, uint64_t>> _ni_connections;
    fz::SafePtr<fz::SafePtr<Float>> _ivpg_to_stif_fi_part;
    fz::SafePtr<fz::SafePtr<Float>> _ivpg_to_mass_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ivpg_to_ptr_in_a;
    fz::SafePtr<fz::SafePtr<Float>> _ispgi_to_damp_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ispgi_to_ptr_in_a;
    fz::SafePtr<Cmplx*> _vpi_to_ptr_in_b;
    fz::SafePtr<fz::SafePtr<Float>> _ispgv_to_forc_fi_part;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _ispgv_to_ptr_in_b;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _pvi_to_ptr_in_a;
    fz::SafePtr<fz::SafePtr<Cmplx*>> _pvi_to_ptr_in_b;
    fz::SafePtr<Cmplx> _a_vals;
    fz::SafePtr<uint64_t> _b_row_idx;
    fz::SafePtr<Cmplx> _b_vals;
    fz::SafePtr<Cmplx> _x;

    uint64_t _fi_count;
    uint64_t _ni_count;
    uint64_t _vei_count;
    uint64_t _isei_count;
    uint64_t _vsei_count;
    uint64_t _psei_count;
    uint64_t _ivpg_count;
    uint64_t _ispgi_count;
    uint64_t _ispgp_count;
    uint64_t _ispgv_count;
    uint64_t _vpi_count;
    uint64_t _pni_count;
    uint64_t _pvi_count;

    std::unique_ptr<indicators::ProgressBar> _progress_bar;

    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        LdltSolver<Cmplx> _solver;
        fz::SafePtr<Cmplx> _b_dense;
        fz::SafePtr<Cmplx> _a_diag;
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        std::optional<Eigen::Map<
            Eigen::SparseMatrix<Cmplx, Eigen::ColMajor, Eigen::Index>
        >> _a_eigen;
        std::array<Eigen::Index, 2UL> _b_col_idx_signed;
        std::optional<Eigen::Map<
            Eigen::SparseMatrix<Cmplx, Eigen::ColMajor, Eigen::Index>
        >> _b_eigen;
        std::optional<Eigen::Map<
            Eigen::Matrix<Cmplx, Eigen::Dynamic, 1UL>
        >> _x_eigen;
        fz::SafePtr<Eigen::Index> _a_row_idx;
        fz::SafePtr<Eigen::Index> _a_col_idx;
        fz::SafePtr<Eigen::Index> _b_row_idx_signed;
        std::unique_ptr<
            Eigen::SparseLU<
                Eigen::SparseMatrix<Cmplx, Eigen::ColMajor, Eigen::Index>,
                Eigen::COLAMDOrdering<Eigen::Index>
            >
        > _solver;
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        fz::SafePtr<Cmplx> _b_dense;
        _MKL_DSS_HANDLE_t _dss_handle;
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        ZMUMPS_STRUC_C _solver;
        fz::SafePtr<MUMPS_INT> _a_row_idx;
        fz::SafePtr<MUMPS_INT> _a_col_idx;
        fz::SafePtr<Cmplx> _b_dense;
    #endif
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
    const char* const hdf5_file_path
);

// alias for simulation type
template<ElementOrder O>
using SimulationFemHelmTet = Simulation<
    NumericalMethod::FEM,
    Equation::HELMHOLTZ,
    ElementShape::TETRAHEDRON,
    O
>;

// macro for explicit instantiation declarations
#define NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3 \
    namespace numav { \
        template class Simulation< \
            NumericalMethod::FEM, \
            Equation::HELMHOLTZ, \
            ElementShape::TETRAHEDRON, \
            ElementOrder::LINEAR \
        >; \
        template class Simulation< \
            NumericalMethod::FEM, \
            Equation::HELMHOLTZ, \
            ElementShape::TETRAHEDRON, \
            ElementOrder::QUADRATIC \
        >; \
    }

} // namespace numav
