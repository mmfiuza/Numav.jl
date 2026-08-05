# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file defines all the subtypes of Simulation.

export Simulation

abstract type Simulation end

@enum ElementShape ElementShape_TETRAHEDRON
@enum ElementOrder ElementOrder_LINEAR ElementOrder_QUADRATIC

function _enis_count(element_order::ElementOrder)
    if element_order === ElementOrder_LINEAR
        return 3
    elseif element_order === ElementOrder_QUADRATIC
        return 4
    end
end

function _eniv_count(element_order::ElementOrder)
    if element_order === ElementOrder_LINEAR
        return 4
    elseif element_order === ElementOrder_QUADRATIC
        return 10
    end
end

@kwdef mutable struct SimulationFemHelmholtz{
    ELEMENT_SHAPE,
    ELEMENT_ORDER,
    ENIS_COUNT,
    ENIV_COUNT,
} <: Simulation
    _fi_to_freq::Vector{Float64} = []
    
    _ni_count::Int = 0
    _sei_count::Int = 0
    _vei_count::Int = 0

    _ni_to_coords::Vector{NTuple{3,Float64}} = []
    _sei_to_ni::Vector{NTuple{ENIS_COUNT,Int}} = []
    _vei_to_ni::Vector{NTuple{ENIV_COUNT,Int}} = []
    _sei_to_espg::Vector{Int} = []
    _vei_to_evpg::Vector{Int} = []

    _cpp_simulation::_cpp_Simulation
end

