# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using Bijections: Bijection

@enum ElementShape ElementShape_TETRAHEDRON
@enum ElementOrder ElementOrder_LINEAR ElementOrder_QUADRATIC

@kwdef mutable struct SimulationFemHelmholtz{
    ELEMENT_SHAPE,
    ELEMENT_ORDER,
    ENIS_COUNT,
    ENIV_COUNT,
} <: Simulation
    _fi_to_freq::Vector{Float64} = []
    
    _fi_count::Int = 0
    _ni_count::Int = 0
    _sei_count::Int = 0
    _vei_count::Int = 0
    _isei_count::Int = 0
    _vsei_count::Int = 0
    _psei_count::Int = 0
    _ivpg_count::Int = 0
    _ispgi_count::Int = 0
    _ispgp_count::Int = 0
    _ispgv_count::Int = 0
    _vpi_count::Int = 0
    _ppi_count::Int = 0
    _pni_count::Int = 0
    _pvi_count::Int = 0

    _ni_to_xyz::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    _sei_to_ni::Matrix{Int} = Matrix{Int}(undef, 0, 0)
    _vei_to_ni::Matrix{Int} = Matrix{Int}(undef, 0, 0)
    _sei_to_espg::Vector{Int} = []
    _vei_to_evpg::Vector{Int} = []

    _vei_to_ivpg::Vector{Int} = []
    _isei_to_ni::Matrix{Int} = Matrix{Int}(undef, 0, 0)
    _isei_to_ispgi::Vector{Int} = []
    _vsei_to_ni::Matrix{Int} = Matrix{Int}(undef, 0, 0)
    _vsei_to_ispgv::Vector{Int} = []
    _psei_to_ni::Matrix{Int} = Matrix{Int}(undef, 0, 0)
    _psei_to_ispgp::Vector{Int} = []

    _existing_evpg::Set{Int} = Set{Int}()
    _existing_espg::Set{Int} = Set{Int}()

    _evpg_ivpg_bimap::Bijection{Int,Int} = Bijection{Int,Int}()
    _espg_ispgi_bimap::Bijection{Int,Int} = Bijection{Int,Int}()
    _espg_ispgv_bimap::Bijection{Int,Int} = Bijection{Int,Int}()
    _espg_ispgp_bimap::Bijection{Int,Int} = Bijection{Int,Int}()

    _ivpg_to_density_func::Vector{Function} = []
    _ivpg_to_soundspeed_func::Vector{Function} = []
    _ispgi_to_impedance_func::Vector{Function} = []
    _ispgv_to_velocity_func::Vector{Function} = []
    _ispgp_to_pressure_func::Vector{Function} = []
    _vpi_to_volvel_func::Vector{Function} = []
    _ppi_to_pressure_func::Vector{Function} = []
    _pvi_to_pressure_func::Vector{Function} = []

    _ivpg_to_density_values::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, 0, 0)
    _ivpg_to_soundspeed_values::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, 0, 0)
    _ispgi_to_impedance_values::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, 0, 0)
    _ispgv_to_velocity_values::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, 0, 0)
    _vpi_to_volvel_values::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, 0, 0)
    _pvi_to_pressure_values::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, 0, 0)

    _vpi_to_ni::Vector{Int} = []
    _ppi_to_ni::Vector{Int} = []
    _pni_to_ni::Vector{Int} = []

    _pvi_to_pni_count::Vector{Int} = []

    _is_freq_defined::Bool = false
    _is_mesh_defined::Bool = false
    _is_any_source_defined::Bool = false
    _did_run::Bool = false

    _hdf5_file_path::String = ""
end

function _enis_count(element_order::ElementOrder)
    if element_order === ElementOrder_LINEAR
        return 3
    elseif element_order === ElementOrder_QUADRATIC
        return 6
    end
end

function _eniv_count(element_order::ElementOrder)
    if element_order === ElementOrder_LINEAR
        return 4
    elseif element_order === ElementOrder_QUADRATIC
        return 10
    end
end

function _enis_count(s::SimulationFemHelmholtz)
    return typeof(s).parameters[3]
end

function _eniv_count(s::SimulationFemHelmholtz)
    return typeof(s).parameters[4]
end

function _is_quadratic(s::SimulationFemHelmholtz)
    return typeof(s).parameters[2] === ElementOrder_QUADRATIC
end
