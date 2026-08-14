# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file defines all the subtypes of Simulation.

export Simulation

abstract type Simulation end

@enum ElementShape ElementShape_TETRAHEDRON
@enum ElementOrder ElementOrder_LINEAR ElementOrder_QUADRATIC

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

    _ni_to_xyz::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    _sei_to_ni::Matrix{Int} = Matrix{Int}(undef, 0, 0)
    _vei_to_ni::Matrix{Int} = Matrix{Int}(undef, 0, 0)
    _sei_to_espg::Vector{Int} = []
    _vei_to_evpg::Vector{Int} = []

    _cpp_simulation::_cpp_Simulation
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
