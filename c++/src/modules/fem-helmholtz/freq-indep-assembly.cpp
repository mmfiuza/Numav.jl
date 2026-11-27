// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "numav/numav.hpp"
#include "modules/fem-helmholtz/fem-helmholtz.hpp"

#include "common/exception.hpp"
#include "common/maths.hpp"
#include "common/hash-functions.hpp"
#include "common/utils.hpp"
#include "modules/fem-helmholtz/shape-functions.hpp"
#include "modules/fem-helmholtz/analytic-integration.hpp"
#include "modules/fem-helmholtz/gauss-quadrature.hpp"

#include <unordered_set>
#include <unordered_map>

namespace numav {

template<typename T>
bool compare_pair(const std::pair<T,T> a, const std::pair<T,T> b) {
    if constexpr
    (GLOBAL_MATRIX_STORAGE_ORDER == MatrixStorageOrder::COL_MAJOR)
    {
        if (a.second != b.second) {
            return a.second < b.second;
        }
        else {
            return a.first < b.first;
        }
    }
    else if constexpr
    (GLOBAL_MATRIX_STORAGE_ORDER == MatrixStorageOrder::ROW_MAJOR)
    {
        if (a.first != b.first) {
            return a.first < b.first;
        }
        else {
            return a.second < b.second;
        }
    }
}

template<typename T>
std::pair<T,T> make_ordered_rowcol_pair(const T a, const T b)
{
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        return std::make_pair(a,b);
    #endif

    if constexpr
    (GLOBAL_MATRIX_TRIANGULAR_TYPE == TriangularMatrixType::UPPER)
    {
        return a<b ? std::make_pair(a,b) : std::make_pair(b,a);
    }
    else if constexpr
    (GLOBAL_MATRIX_TRIANGULAR_TYPE == TriangularMatrixType::LOWER)
    {
        return a<b ? std::make_pair(b,a) : std::make_pair(a,b);
    }
}

template <ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_allocate_a()
{
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        constexpr std::array<
            std::array<uint64_t, 2UL>,
            PERMUTATION_REP_SIZE<ENIV_COUNT<S,O>, 2UL>
        > ENI_PAIRS = PERMUTATION_REP<ENIV_COUNT<S,O>>;
    #else
        constexpr std::array<
            std::array<uint64_t, 2UL>,
            COMBINATION_REP_SIZE<ENIV_COUNT<S,O>, 2UL>
        > ENI_PAIRS = COMBINATION_REP<ENIV_COUNT<S,O>>;
    #endif
    std::unordered_set<std::pair<uint64_t,uint64_t>> existing_pairs;
    for (uint64_t vei = 0UL; vei != _vei_count; ++vei) {
        for (const auto& eniv : ENI_PAIRS) {
            #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
                if (eniv[0UL] == eniv[1UL]) { continue; }
            #endif
            existing_pairs.insert(
                make_ordered_rowcol_pair(
                    _vei_to_ni[vei][eniv[0UL]], _vei_to_ni[vei][eniv[1UL]]
                )
            );
        }
    }
    _ni_connections = fz::SafePtr<std::pair<uint64_t,uint64_t>>(
        existing_pairs.size()
    );
    std::copy(
        existing_pairs.begin(),
        existing_pairs.end(),
        _ni_connections.begin()
    );
    std::sort(
        _ni_connections.begin(),
        _ni_connections.end(),
        compare_pair<uint64_t>
    );
    _a_vals = fz::SafePtr<Cmplx>(_ni_connections.size());
    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        _a_diag = fz::SafePtr<Cmplx>(_ni_count);
    #endif
}

template <ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_allocate_b()
{
    std::unordered_set<uint64_t> existing_source_nodes;
    for (uint64_t vsei = 0UL; vsei != _vsei_count; ++vsei) {
        for (uint64_t enis = 0UL; enis != ENIS_COUNT<S,O>; ++enis) {
            const uint64_t ni = _vsei_to_ni[vsei][enis];
            existing_source_nodes.insert(ni);
        }
    }
    for (uint64_t vpi = 0UL; vpi != _vpi_count; ++vpi) {
        existing_source_nodes.insert(_vpi_to_ni[vpi]);
    }
    for (uint64_t pni = 0UL; pni != _pni_count; ++pni) {
        existing_source_nodes.insert(_pni_to_ni[pni]);
    }
    _b_row_idx = fz::SafePtr<uint64_t>(existing_source_nodes.size());
    std::copy(
        existing_source_nodes.begin(),
        existing_source_nodes.end(),
        _b_row_idx.begin()
    );
    std::sort(_b_row_idx.begin(), _b_row_idx.end());
    _b_vals = fz::SafePtr<Cmplx>(_b_row_idx.size());
}

template<ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_assemble_fi_part_for_vol_elements()
{
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        constexpr std::array<
            std::array<uint64_t, 2UL>,
            PERMUTATION_REP_SIZE<ENIV_COUNT<S,O>, 2UL>
        > ENI_PAIRS = PERMUTATION_REP<ENIV_COUNT<S,O>>;
    #else
        constexpr std::array<
            std::array<uint64_t, 2UL>,
            COMBINATION_REP_SIZE<ENIV_COUNT<S,O>, 2UL>
        > ENI_PAIRS = COMBINATION_REP<ENIV_COUNT<S,O>>;
    #endif

    // count the fipi for each ivr
    fz::SafePtr<std::unordered_map<
        std::pair<uint64_t,uint64_t>, uint64_t
    >> ivr_to_map_to_fipi(_ivr_count);
    for (uint64_t vei = 0UL; vei != _vei_count; ++vei)
    {
        const uint64_t ivr = _vei_to_ivr[vei];
        for (const auto& eni : ENI_PAIRS) {
            const std::pair<uint64_t,uint64_t> pair = make_ordered_rowcol_pair(
                _vei_to_ni[vei][eni[0UL]], _vei_to_ni[vei][eni[1UL]]
            );
            if (!(ivr_to_map_to_fipi[ivr].count(pair) > 0)) {
                const uint64_t fipi = ivr_to_map_to_fipi[ivr].size();
                ivr_to_map_to_fipi[ivr].insert({pair, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _ivr_to_stif_fi_part = fz::SafePtr<fz::SafePtr<Float>>(_ivr_count);
    _ivr_to_mass_fi_part = fz::SafePtr<fz::SafePtr<Float>>(_ivr_count);
    _ivr_to_ptr_in_a = fz::SafePtr<fz::SafePtr<Cmplx*>>(_ivr_count);
    for (uint64_t ivr = 0UL; ivr != _ivr_count; ++ivr)
    {
        const uint64_t size = ivr_to_map_to_fipi[ivr].size();
        _ivr_to_stif_fi_part[ivr] = fz::SafePtr<Float>(size);
        _ivr_to_stif_fi_part[ivr].fill(0_F);
        _ivr_to_mass_fi_part[ivr] = fz::SafePtr<Float>(size);
        _ivr_to_mass_fi_part[ivr].fill(0_F);
        _ivr_to_ptr_in_a[ivr] = fz::SafePtr<Cmplx*>(size);
    }

    // assemble elemental stiffness and mass matrices
    std::array<uint64_t, ENI_PAIRS.size()> fipi_vol;
    for (uint64_t vei = 0UL; vei != _vei_count; ++vei)
    {
        const uint64_t ivr = _vei_to_ivr[vei];
        
        // Create _ivr_to_ptr_in_a
        for (uint64_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci)
        {
            const std::pair<uint64_t,uint64_t> pair = make_ordered_rowcol_pair(
                _vei_to_ni[vei][ENI_PAIRS[nci][0UL]],
                _vei_to_ni[vei][ENI_PAIRS[nci][1UL]]
            );
            fipi_vol[nci] = ivr_to_map_to_fipi[ivr].at(pair);

            #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
                if(pair.first == pair.second) {
                    const uint64_t ni = pair.first;
                    _ivr_to_ptr_in_a[ivr][fipi_vol[nci]] =
                        _a_diag.data() + ni;
                    continue;
                }
            #endif

            const std::pair<uint64_t,uint64_t>* const pair_ptr =
            std::lower_bound(
                _ni_connections.begin(),
                _ni_connections.end(),
                pair,
                compare_pair<uint64_t>
            );
            const uint64_t ptrdiff = pair_ptr - _ni_connections.begin();
            _ivr_to_ptr_in_a[ivr][fipi_vol[nci]] = _a_vals.begin() + ptrdiff;
        }
        
        // coordinates matrix
        Eigen::Matrix<Float,DIM,ENIV_COUNT<S,O>> coords_matrix;
        for (uint64_t eni = 0UL; eni != ENIV_COUNT<S,O>; ++eni) {
            const uint64_t ni = _vei_to_ni[vei][eni];
            for (uint64_t di = 0UL; di != DIM; ++di) {
                coords_matrix(di,eni) = _ni_to_xyz[ni][di];
            }
        }
        #if \
        NUMAV_STIF_INTEGRATION_METHOD == NUMAV_ANALYTIC || \
        NUMAV_MASS_INTEGRATION_METHOD == NUMAV_ANALYTIC
            const Float tet_volume = get_tetrahedron_volume(
                *reinterpret_cast<std::array<std::array<Float,3UL>,4UL>*>(
                    coords_matrix.data()
                )
            );
        #endif

        // stiffness matrix
        #if NUMAV_STIF_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            constexpr std::array<std::array<Float,DIM>,NGP_STIF<S,O>>
                GAUSS_POINTS_STIF = GAUSS_POINTS_VOL<S,NGP_STIF<S,O>>;
            for (uint64_t gpi = 0UL; gpi != NGP_STIF<S,O>; ++gpi)
            {
                const Eigen::Matrix<Float,ENIV_COUNT<S,O>,DIM> nabla_n =
                    shape_func_vol_gradient<S,O>(
                        GAUSS_POINTS_STIF[gpi][0UL],
                        GAUSS_POINTS_STIF[gpi][1UL],
                        GAUSS_POINTS_STIF[gpi][2UL]
                    );

                const Eigen::Matrix<Float,DIM,DIM> jac_matrix = 
                    coords_matrix * nabla_n;
                
                const Float detj = std::abs(jac_matrix.determinant());

                const Eigen::Matrix<Float,DIM,DIM> inv_jac =
                    jac_matrix.inverse();
                
                const Eigen::Matrix<Float,ENIV_COUNT<S,O>,DIM> b_matrix =
                    nabla_n * inv_jac;
                
                Eigen::Matrix<Float,ENIV_COUNT<S,O>,ENIV_COUNT<S,O>> bbt =
                    b_matrix * b_matrix.transpose();
                
                auto& bbt_detj_w = bbt;
                bbt_detj_w = bbt*detj*GAUSS_WEIGHTS_VOL<S,NGP_STIF<S,O>>[gpi];

                for (uint64_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                    _ivr_to_stif_fi_part[ivr][fipi_vol[nci]] +=
                        bbt_detj_w(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
                }
            }
        #elif NUMAV_STIF_INTEGRATION_METHOD == NUMAV_ANALYTIC
            Eigen::Matrix<Float, ENIV_COUNT<S,O>, ENIV_COUNT<S,O>> stif_fi_part = 
                get_stif_matrix_const_part<S,O>(coords_matrix);

            stif_fi_part /= tet_volume;

            for (uint64_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                _ivr_to_stif_fi_part[ivr][fipi_vol[nci]] +=
                    stif_fi_part(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
            }
        #endif

        // mass matrix
        #if NUMAV_MASS_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            constexpr std::array<std::array<Float,DIM>,NGP_MASS<S,O>>
                GAUSS_POINTS_MASS = GAUSS_POINTS_VOL<S,NGP_MASS<S,O>>;
            for (uint64_t gpi = 0UL; gpi != NGP_MASS<S,O>; ++gpi)
            {
                const Eigen::Matrix<Float,ENIV_COUNT<S,O>,DIM> nabla_n =
                    shape_func_vol_gradient<S,O>(
                        GAUSS_POINTS_MASS[gpi][0UL],
                        GAUSS_POINTS_MASS[gpi][1UL],
                        GAUSS_POINTS_MASS[gpi][2UL]
                    );

                const Eigen::Matrix<Float,DIM,DIM> jac_matrix = 
                    coords_matrix * nabla_n;
                
                const Float detj = std::abs(jac_matrix.determinant());
                    
                const Eigen::Matrix<Float,ENIV_COUNT<S,O>,1UL> n =
                    shape_func_vol<S,O>(
                        GAUSS_POINTS_MASS[gpi][0UL],
                        GAUSS_POINTS_MASS[gpi][1UL],
                        GAUSS_POINTS_MASS[gpi][2UL]
                    );
                    
                Eigen::Matrix<Float,ENIV_COUNT<S,O>,ENIV_COUNT<S,O>> nnt =
                    n * n.transpose();
                    
                auto& nnt_detj_w = nnt; 
                nnt_detj_w = nnt*detj*GAUSS_WEIGHTS_VOL<S,NGP_MASS<S,O>>[gpi];
                
                for (uint64_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                    _ivr_to_mass_fi_part[ivr][fipi_vol[nci]] -=
                        nnt_detj_w(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
                }
            }
        #elif NUMAV_MASS_INTEGRATION_METHOD == NUMAV_ANALYTIC
            Eigen::Matrix<Float,ENIV_COUNT<S,O>,ENIV_COUNT<S,O>> mass_fi_part =
                MASS_MATRIX_CONST_PART<S,O>;
            
            mass_fi_part *= tet_volume;
            
            for (uint64_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                _ivr_to_mass_fi_part[ivr][fipi_vol[nci]] -=
                    mass_fi_part(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
            }
        #endif
    }
    ivr_to_map_to_fipi.free();
}

template<ElementShape S, ElementOrder O>
std::array<std::array<Float,2UL>,ENIS_COUNT<S,O>> project_triangle_to_2d(
    const std::array<std::array<Float,DIM>,ENIS_COUNT<S,O>> vertices_3d
) {
    std::array<Eigen::Vector3d, ENIS_COUNT<S,O>> point_minus_origin;
    for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
        point_minus_origin[eni] = Eigen::Vector3d(
            vertices_3d[eni][0UL] - vertices_3d[0UL][0UL],
            vertices_3d[eni][1UL] - vertices_3d[0UL][1UL],
            vertices_3d[eni][2UL] - vertices_3d[0UL][2UL]
        );
    }
    const Eigen::Vector3d& u = point_minus_origin[1UL];
    const Eigen::Vector3d& v = point_minus_origin[2UL];
    const Eigen::Vector3d n = u.cross(v);
    const Eigen::Vector3d x = u / u.norm();
    Eigen::Vector3d y = n.cross(u);
    y /= y.norm();

    std::array<std::array<Float,2UL>,ENIS_COUNT<S,O>> vertices_2d;
    for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
        vertices_2d[eni] = std::array<Float,2UL>({
            point_minus_origin[eni].dot(x), point_minus_origin[eni].dot(y)
        });
    }
    return vertices_2d;
}

template<ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_assemble_fi_part_for_sfc_impedance()
{
    #if NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        constexpr std::array<
            std::array<uint64_t, 2UL>,
            PERMUTATION_REP_SIZE<ENIS_COUNT<S,O>, 2UL>
        > ENI_PAIRS = PERMUTATION_REP<ENIS_COUNT<S,O>>;
    #else
        constexpr std::array<
            std::array<uint64_t, 2UL>,
            COMBINATION_REP_SIZE<ENIS_COUNT<S,O>, 2UL>
        > ENI_PAIRS = COMBINATION_REP<ENIS_COUNT<S,O>>;
    #endif

    // count the fipi for each isri
    fz::SafePtr<std::unordered_map<
        std::pair<uint64_t,uint64_t>, uint64_t
    >> isri_to_map_to_fipi(_isri_count);
    for (uint64_t isei = 0UL; isei != _isei_count; ++isei)
    {
        const uint64_t isri = _isei_to_isri[isei];
        for (const auto& eni : ENI_PAIRS) {
            const std::pair<uint64_t,uint64_t> pair = make_ordered_rowcol_pair(
                _isei_to_ni[isei][eni[0UL]], _isei_to_ni[isei][eni[1UL]]
            );
            if (!(isri_to_map_to_fipi[isri].count(pair) > 0)) {
                const uint64_t fipi = isri_to_map_to_fipi[isri].size();
                isri_to_map_to_fipi[isri].insert({pair, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _isri_to_damp_fi_part = fz::SafePtr<fz::SafePtr<Float>>(_isri_count);
    _isri_to_ptr_in_a = fz::SafePtr<fz::SafePtr<Cmplx*>>(_isri_count);
    for (uint64_t isri = 0UL; isri != _isri_count; ++isri)
    {
        const uint64_t size = isri_to_map_to_fipi[isri].size();
        _isri_to_damp_fi_part[isri] = fz::SafePtr<Float>(size);
        _isri_to_damp_fi_part[isri].fill(0_F);
        _isri_to_ptr_in_a[isri] = fz::SafePtr<Cmplx*>(size);
    }

    // assemble the elemental damping matrices
    std::array<uint64_t, ENI_PAIRS.size()> fipi_damp;
    for (uint64_t isei = 0UL; isei != _isei_count; ++isei)
    {
        const uint64_t isri = _isei_to_isri[isei];
        
        // Create _isri_to_ptr_in_a
        for (uint64_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci)
        {
            const std::pair<uint64_t,uint64_t> pair = make_ordered_rowcol_pair(
                _isei_to_ni[isei][ENI_PAIRS[nci][0UL]],
                _isei_to_ni[isei][ENI_PAIRS[nci][1UL]]
            );
            fipi_damp[nci] = isri_to_map_to_fipi[isri].at(pair);
            
            #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
                if(pair.first == pair.second) {
                    const uint64_t ni = pair.first;
                    _isri_to_ptr_in_a[isri][fipi_damp[nci]] =
                        _a_diag.data() + ni;
                    continue;
                }
            #endif

            const std::pair<uint64_t,uint64_t>* const pair_ptr =
            std::lower_bound(
                _ni_connections.begin(),
                _ni_connections.end(),
                pair,
                compare_pair<uint64_t>
            );
            const uint64_t ptrdiff = pair_ptr - _ni_connections.begin();
            _isri_to_ptr_in_a[isri][fipi_damp[nci]] =
                _a_vals.begin() + ptrdiff;
        }

        // damping matrix
        #if NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            std::array<std::array<Float,DIM>, ENIS_COUNT<S,O>> triangle_3d;
            for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
                const uint64_t ni = _isei_to_ni[isei][eni];
                triangle_3d[eni] = _ni_to_xyz[ni];
            }
            std::array<std::array<Float,2UL>,ENIS_COUNT<S,O>> triangle_2d =
                project_triangle_to_2d<S,O>(triangle_3d);

            Eigen::Matrix<Float, 2UL, ENIS_COUNT<S,O>> coords_matrix;
            for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
                coords_matrix(0UL,eni) = triangle_2d[eni][0UL];
                coords_matrix(1UL,eni) = triangle_2d[eni][1UL];
            }

            constexpr std::array<std::array<Float,2UL>,NGP_DAMP<S,O>>
                GAUSS_POINTS_DAMP = GAUSS_POINTS_SFC<S,NGP_DAMP<S,O>>;
            for (uint64_t gpi = 0UL; gpi != NGP_DAMP<S,O>; ++gpi)
            {
                const Eigen::Matrix<Float,ENIS_COUNT<S,O>,2UL> nabla_n =
                    shape_func_sfc_gradient<S,O>(
                        GAUSS_POINTS_DAMP[gpi][0UL],
                        GAUSS_POINTS_DAMP[gpi][1UL]
                    );
                
                const Eigen::Matrix<Float,2UL,2UL> jac_matrix =
                    coords_matrix * nabla_n;
                
                const Float detj = std::abs(jac_matrix.determinant());
                
                const Eigen::Matrix<Float,ENIS_COUNT<S,O>,1UL> n =
                    shape_func_sfc<S,O>(
                        GAUSS_POINTS_DAMP[gpi][0UL],
                        GAUSS_POINTS_DAMP[gpi][1UL]
                    );
                
                Eigen::Matrix<Float,ENIS_COUNT<S,O>,ENIS_COUNT<S,O>> nnt =
                    n * n.transpose();
                                
                auto& nnt_detj_w = nnt;
                nnt_detj_w = nnt*detj*GAUSS_WEIGHTS_SFC<S,NGP_DAMP<S,O>>[gpi];
                
                for (uint64_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                    _isri_to_damp_fi_part[isri][fipi_damp[nci]] += 
                        nnt_detj_w(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
                }
            }
        #elif NUMAV_DAMP_INTEGRATION_METHOD == NUMAV_ANALYTIC
            std::array<std::array<Float,DIM>,3UL> vertex_coords;
            for (uint64_t eni = 0UL; eni != 3UL; ++eni) {
                const uint64_t ni = _isei_to_ni[isei][eni];
                vertex_coords[eni] = _ni_to_xyz[ni];
            }
            const Float triangle_area = get_triangle_area(vertex_coords);

            Eigen::Matrix<Float,ENIS_COUNT<S,O>,ENIS_COUNT<S,O>> damp_fi_part =
                DAMP_MATRIX_CONST_PART<S,O>;
            
            damp_fi_part *= triangle_area;
            
            for (uint64_t nci = 0UL; nci != ENI_PAIRS.size(); ++nci) {
                _isri_to_damp_fi_part[isri][fipi_damp[nci]] +=
                    damp_fi_part(ENI_PAIRS[nci][0UL], ENI_PAIRS[nci][1UL]);
            }
        #else
            static_assert(false, "Invalid NUMAV_DAMP_INTEGRATION_METHOD.");
        #endif
        }
    isri_to_map_to_fipi.free();
}

template<ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_assemble_fi_part_for_sfc_velocity()
{
    // count the fipi for each isrv
    fz::SafePtr<std::unordered_map<uint64_t,uint64_t>> isrv_to_map_to_fipi(
        _isrv_count
    );
    for (uint64_t vsei = 0UL; vsei != _vsei_count; ++vsei) {
        const uint64_t isrv = _vsei_to_isrv[vsei];
        for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
            const uint64_t ni = _vsei_to_ni[vsei][eni];
            if (!(isrv_to_map_to_fipi[isrv].count(ni) > 0)) {
                const uint64_t fipi = isrv_to_map_to_fipi[isrv].size();
                isrv_to_map_to_fipi[isrv].insert({ni, fipi});
            }
        }
    }

    // allocate memory in the safe pointers
    _isrv_to_forc_fi_part = fz::SafePtr<fz::SafePtr<Float>>(_isrv_count);
    _isrv_to_ptr_in_b = fz::SafePtr<fz::SafePtr<Cmplx*>>(_isrv_count);
    for (uint64_t isrv = 0UL; isrv != _isrv_count; ++isrv)
    {
        const uint64_t size = isrv_to_map_to_fipi[isrv].size();
        _isrv_to_forc_fi_part[isrv] = fz::SafePtr<Float>(size);
        _isrv_to_forc_fi_part[isrv].fill(0_F);
        _isrv_to_ptr_in_b[isrv] = fz::SafePtr<Cmplx*>(size);
    }

    // assemble elemental force vectors
    std::array<uint64_t, ENIS_COUNT<S,O>> fipi_forc;
    for (uint64_t vsei = 0UL; vsei != _vsei_count; ++vsei)
    {
        const uint64_t isrv = _vsei_to_isrv[vsei];

        // Create _isrv_to_ptr_in_b
        for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni)
        {
            const uint64_t ni = _vsei_to_ni[vsei][eni];
            fipi_forc[eni] = isrv_to_map_to_fipi[isrv].at(ni);
            
            const uint64_t* const ni_ptr = std::lower_bound(
                _b_row_idx.begin(), _b_row_idx.end(), ni
            );
            const uint64_t ptrdiff = ni_ptr - _b_row_idx.begin();
            _isrv_to_ptr_in_b[isrv][fipi_forc[eni]] =
                _b_vals.begin() + ptrdiff;
        }

        // elemental force vector
        #if NUMAV_FORC_INTEGRATION_METHOD == NUMAV_GAUSS_QUADRATURE
            std::array<std::array<Float,DIM>, ENIS_COUNT<S,O>> triangle_3d;
            for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
                const uint64_t ni = _vsei_to_ni[vsei][eni];
                triangle_3d[eni] = _ni_to_xyz[ni];
            }
            std::array<std::array<Float,2UL>,ENIS_COUNT<S,O>> triangle_2d =
                project_triangle_to_2d<S,O>(triangle_3d);

            Eigen::Matrix<Float, 2UL, ENIS_COUNT<S,O>> coords_matrix;
            for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
                coords_matrix(0UL,eni) = triangle_2d[eni][0UL];
                coords_matrix(1UL,eni) = triangle_2d[eni][1UL];
            }

            constexpr std::array<std::array<Float,2UL>,NGP_FORC<S,O>>
                GAUSS_POINTS_FORC = GAUSS_POINTS_SFC<S,NGP_FORC<S,O>>;
            for (uint64_t gpi = 0UL; gpi != NGP_FORC<S,O>; ++gpi)
            {
                const Eigen::Matrix<Float,ENIS_COUNT<S,O>, 2UL> nabla_n =
                    shape_func_sfc_gradient<S,O>(
                        GAUSS_POINTS_FORC[gpi][0UL],
                        GAUSS_POINTS_FORC[gpi][1UL]
                    );
                
                const Eigen::Matrix<Float,2UL,2UL> jac_matrix =
                    coords_matrix * nabla_n;
                
                const Float detj = std::abs(jac_matrix.determinant());

                Eigen::Matrix<Float,ENIS_COUNT<S,O>,1UL> n =
                    shape_func_sfc<S,O>(
                        GAUSS_POINTS_FORC[gpi][0UL],
                        GAUSS_POINTS_FORC[gpi][1UL]
                    );

                auto& n_detj_w = n;
                n_detj_w = n * detj * GAUSS_WEIGHTS_SFC<S,NGP_FORC<S,O>>[gpi];
                
                for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
                    _isrv_to_forc_fi_part[isrv][fipi_forc[eni]] +=
                        n_detj_w(eni);
                }
            }
        #elif NUMAV_FORC_INTEGRATION_METHOD == NUMAV_ANALYTIC
            std::array<std::array<Float,DIM>,3UL> vertex_coords;
            for (uint64_t eni = 0UL; eni != 3UL; ++eni) {
                const uint64_t ni = _vsei_to_ni[vsei][eni];
                vertex_coords[eni] = _ni_to_xyz[ni];
            }
            const Float triangle_area = get_triangle_area(vertex_coords);

            Eigen::Matrix<Float,ENIS_COUNT<S,O>,1UL> forc_fi_part =
                FORC_VECTOR_CONST_PART<S,O>;
            
            forc_fi_part *= triangle_area;
            
            for (uint64_t eni = 0UL; eni != ENIS_COUNT<S,O>; ++eni) {
                _isrv_to_forc_fi_part[isrv][fipi_forc[eni]] +=
                    forc_fi_part(eni);
            }
        #else
            static_assert(false, "Invalid NUMAV_FORC_INTEGRATION_METHOD.");
        #endif
    }
    isrv_to_map_to_fipi.free();
}

template<ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_assemble_fi_part_for_point_velocity()
{
    _vpi_to_ptr_in_b = fz::SafePtr<Cmplx*>(_vpi_count);
    for (uint64_t vpi = 0UL; vpi != _vpi_count; ++vpi)
    {
        const uint64_t ni = _vpi_to_ni[vpi];
        const uint64_t* const ni_ptr = std::lower_bound(
            _b_row_idx.begin(), _b_row_idx.end(), ni
        );
        const uint64_t ptrdiff = ni_ptr - _b_row_idx.begin();
        _vpi_to_ptr_in_b[vpi] = _b_vals.begin() + ptrdiff;
    }
}

template<ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_assemble_fi_part_for_pressure()
{
    // define _pvi_to_ptr_in_a and _pvi_to_ptr_in_b
    _pvi_to_ptr_in_a = fz::SafePtr<fz::SafePtr<Cmplx*>>(_pvi_count);
    _pvi_to_ptr_in_b = fz::SafePtr<fz::SafePtr<Cmplx*>>(_pvi_count);
    uint64_t pni = 0UL;
    for (uint64_t pvi = 0UL; pvi != _pvi_count; ++pvi) {
        const uint64_t pni_count = _pvi_to_pni_count[pvi];
        _pvi_to_ptr_in_a[pvi] = fz::SafePtr<Cmplx*>(pni_count);
        _pvi_to_ptr_in_b[pvi] = fz::SafePtr<Cmplx*>(pni_count);
        for (uint64_t fipi = 0UL; fipi != pni_count; ++fipi) {
            const uint64_t ni = _pni_to_ni[pni];
            ++pni;

            // _pvi_to_ptr_in_a
            #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
                _pvi_to_ptr_in_a[pvi][fipi] = _a_diag.data() + ni;
            #else
                const std::pair<uint64_t,uint64_t> pair = 
                    std::make_pair(ni, ni);
                const std::pair<uint64_t,uint64_t>* const pair_ptr =
                    std::lower_bound(
                        _ni_connections.begin(),
                        _ni_connections.end(),
                        pair,
                        compare_pair<uint64_t>
                    );
                const uint64_t ptrdiff_a = pair_ptr - _ni_connections.begin();
                _pvi_to_ptr_in_a[pvi][fipi] = _a_vals.begin() + ptrdiff_a;
            #endif

            // _pvi_to_ptr_in_b
            const uint64_t* const ni_ptr = std::lower_bound(
                _b_row_idx.begin(), _b_row_idx.end(), ni
            );
            const uint64_t ptrdiff_b = ni_ptr - _b_row_idx.begin();
            _pvi_to_ptr_in_b[pvi][fipi] = _b_vals.begin() + ptrdiff_b;
        }
    }
}

template<ElementShape S, ElementOrder O>
void SimulationFemHelmholtz<S,O>::_assemble_freq_independent_parts()
{   
    _allocate_a();
    _allocate_b();
    #if NUMAV_SYSTEM_SOLVER == NUMAV_INTERNAL
        _define_sparsity_pattern_using_internal_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_EIGEN
        _define_sparsity_pattern_using_eigen_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_ONEMKL
        _define_sparsity_pattern_using_onemkl_solver();
    #elif NUMAV_SYSTEM_SOLVER == NUMAV_MUMPS
        _define_sparsity_pattern_using_mumps_solver();
    #else
        static_assert(false, "Invalid NUMAV_SYSTEM_SOLVER.");
    #endif
    _assemble_fi_part_for_vol_elements();
    _assemble_fi_part_for_sfc_impedance();
    _assemble_fi_part_for_point_velocity();
    _assemble_fi_part_for_sfc_velocity();
    _assemble_fi_part_for_pressure();
}

} // namespace numav

NUMAV_INSTANTIATE_SIM_AC_FEM_FREQ_D3
