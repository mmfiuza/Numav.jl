# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    add_volume_material!,
    add_surface_material!,
    add_sound_source!,
    plot_pressure_field

using HDF5
using Plots
using Statistics: mean
using GLMakie
using LinearAlgebra 

@kwdef mutable struct SimulationFemHelmholtz{
    S<:ElementShape,
    O<:ElementOrder,
} <: Simulation

    numerical_method::FemNumericalMethod = FemNumericalMethod()
    equation::HelmholtzEquation = HelmholtzEquation()
    element_shape::S = S()
    element_order::O = O()

    _cpp_simulation::_cpp_Simulation
    _fi_to_freq::Vector{Float64} = []
end

function _check_if_simulation_has_fem_helmholtz(s::Simulation)
    if (
        !hasproperty(s, :numerical_method) ||
        !hasproperty(s, :equation) ||
        !(
            (s.numerical_method isa FemNumericalMethod) &&
            (s.equation isa HelmholtzEquation)
        )
    )
        _throw_simulation_not_applicable()
    end
    return
end

"""
Assigns acoustic material properties to a volumetric region of the mesh, identified by its physical group tag.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| **Keyword arguments** | | |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID (or vector of IDs) from the mesh |
| `density` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Density in kg/m³ |
| `speed_of_sound` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Speed of sound in m/s |

---
# Examples

> ```julia
> rho(f) = 1.20 # air density in kg/m³
> c(f) = 343 # speed of sound in m/s
> add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)
> ```

> Physical quantities can by complex to model sound absorption:
> ```julia
> rho(f) = 1.20 + 0.001*f
> c(f) = 340 + 0.1*sqrt(f)
> add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)
> ```

> Multiple physical groups can be assigned at once:
> ```julia
> add_volume_material!(s, physical_group=[4,6], density=1.2, speed_of_sound=343)
> ```
"""
function add_volume_material!( 
    simulation::SimulationFemHelmholtz;
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    density::Fdpq,
    speed_of_sound::Fdpq
)
    _check_if_simulation_has_fem_helmholtz(simulation)
    if physical_group isa AbstractVector
        for pg in physical_group
            add_volume_material!(
                simulation,
                physical_group = pg,
                density = density,
                speed_of_sound = speed_of_sound
            )
        end
        return
    end
    density = _fdpq_to_function(density)
    speed_of_sound = _fdpq_to_function(speed_of_sound)
    _cpp_add_volume_material!(
        simulation._cpp_simulation,
        UInt64(physical_group),
        _cmplx_split_and_store(density)...,
        _cmplx_split_and_store(speed_of_sound)...
    )
end

"""
Assigns a specific acoustic impedance boundary condition to a surface in the mesh. This is used to model absorbers, reflecting surfaces or other boundary treatments.

The specific acoustic impedance is the ratio of complex amplitude of acoustic pressure to complex amplitude of normal particle velocity (in Pa·s/m).

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| **Keyword arguments** | | |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID of the boundary surface |
| `specific_acoustic_impedance` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | specific surface acoustic impedance (Pa·s/m) |

---
# Examples

> ```julia
> Z(f) = 1f + 2im
> add_surface_material!(s, physical_group=4, specific_acoustic_impedance=Z)
> ```

> Multiple physical groups can be assigned at once:
> ```julia
> add_surface_material!(s, physical_group=[3,1], specific_acoustic_impedance=1.0)
> ```
"""
function add_surface_material!( 
    simulation::SimulationFemHelmholtz;
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    specific_acoustic_impedance::Fdpq
)
    _check_if_simulation_has_fem_helmholtz(simulation)
    if physical_group isa AbstractVector
        for pg in physical_group
            add_surface_material!(
                simulation,
                physical_group = pg,
                specific_acoustic_impedance = specific_acoustic_impedance
            )
        end
        return
    end
    specific_acoustic_impedance = _fdpq_to_function(specific_acoustic_impedance)
    _cpp_add_surface_material!(
        simulation._cpp_simulation,
        UInt64(physical_group),
        _cpp_PhysicalQuantity_impedance,
        _cmplx_split_and_store(specific_acoustic_impedance)...
    )
end

"""
Sources can be applied either at a specific point in space (via `coordinates`) or over a surface region (via `physical_group`). Three excitation types are supported: `volume_velocity`, `particle velocity`, and `pressure`.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| **Keyword arguments** | | |
| `coordinates` | `Vector{Real}`, `Vector{Vector{Real}}` | `[x, y, z]` location of a point source in m |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID of a surface or volume region |
| `volume_velocity` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Volume velocity in m³/s |
| `particle_velocity` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Normal particle velocity in m/s |
| `pressure` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Acoustic pressure in Pa |

---
# Examples

> Volume velocity source (monopole):
> ```julia
> Q(f) = 10/f # Volume velocity in m³/s as a function of frequency
> add_sound_source!(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)
> ```
> Suitable for modeling punctual omnidirectional sources.

> Particle velocity source (vibrating surface):
> ```julia
> U(f) = 15/f # Particle velocity in m/s
> add_sound_source!(s, physical_group=2, particle_velocity=U)
> ```
> Prescribes the normal component of particle velocity on all surfaces of a physical group. Useful for modeling vibrating panels or pistons.

> Pressure source:
> ```julia
> P(f) = 2f # Pressure in Pa as a function of frequency
> 
> # At a point in space
> add_sound_source!(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
> 
> # On a mesh surface
> add_sound_source!(s, physical_group=3, pressure=P)
> ```
> Prescribes acoustic pressure, either at a point or on a surface.

> Multiple points or physical groups can be assigned at once:
> ```julia
> p1 = [1.0, 3.0, 2.0]
> p2 = [3.0, 1.0, 1.0]
> add_sound_source!(s, coordinates=[p1,p2], volume_velocity=0.01)
> add_sound_source!(s, physical_group=[3,5,9,2], particle_velocity=0.01)
> ```

!!! note
    Each call to `add_sound_source!` should specify either `coordinates` or `physical_group` (not both), and exactly **one** excitation type keyword.
"""
function add_sound_source!( 
    simulation::SimulationFemHelmholtz;
    coordinates::Union{
        AbstractVector{<:Real},
        AbstractVector{<:AbstractVector{<:Real}},
        Nothing
    } = nothing,
    physical_group::Union{
        Integer,
        AbstractVector{<:Integer},
        Nothing
    } = nothing,
    volume_velocity::Union{Fdpq, Nothing} = nothing,
    particle_velocity::Union{Fdpq, Nothing} = nothing,
    pressure::Union{Fdpq, Nothing} = nothing
)
    _check_if_simulation_has_fem_helmholtz(simulation)
    if isnothing(coordinates) && isnothing(physical_group)
        throw(ArgumentError(
            "`coordinates` and `physical_group` not defined"
        ))
    end
    if !isnothing(coordinates) && !isnothing(physical_group)
        throw(ArgumentError(
            "`coordinates` and `physical_group` defined simultaneously"
        ))
    end
    fdpq_count = count(
        !isnothing, (volume_velocity, particle_velocity, pressure)
    )
    if fdpq_count == 0
        throw(ArgumentError(
            "`volume_velocity`, `particle_velocity` and `pressure` not defined"
        ))
    end
    if fdpq_count > 1
        throw(ArgumentError(
            "`volume_velocity`, `particle_velocity` or `pressure` defined " *
            "simultaneously"
        ))
    end
    if coordinates isa AbstractVector{<:Real} && length(coordinates) != 3
        throw(ArgumentError(
            "x,y,z coordinates does not have 3 components"
        ))
    end
    if coordinates isa AbstractVector{<:AbstractVector}
        for c in coordinates
            add_sound_source!(
                simulation,
                coordinates = c,
                physical_group = physical_group,
                volume_velocity = volume_velocity,
                particle_velocity = particle_velocity,
                pressure = pressure
            )
        end
        return
    end
    if physical_group isa AbstractVector
        for pg in physical_group
            add_sound_source!(
                simulation,
                coordinates = coordinates,
                physical_group = pg,
                volume_velocity = volume_velocity,
                particle_velocity = particle_velocity,
                pressure = pressure
            )
        end
        return
    end

    # Check if volume_velocity, particle_velocity or pressure was given
    fdpqv = Ref{Fdpq}()
    if !isnothing(volume_velocity)
        pq_type = _cpp_PhysicalQuantity_volume_velocity
        fdpqv[] = volume_velocity
    elseif !isnothing(particle_velocity)
        pq_type = _cpp_PhysicalQuantity_particle_velocity
        fdpqv[] = particle_velocity
    elseif !isnothing(pressure)
        pq_type = _cpp_PhysicalQuantity_pressure
        fdpqv[] = pressure
    end

    fdpqv[] = _fdpq_to_function(fdpqv[])
    source_args =
    if !isnothing(coordinates)
        (_cpp_SourceType_point, Float64.(coordinates))
    elseif !isnothing(physical_group)
        (_cpp_SourceType_surface, UInt64(physical_group))
    end

    _cpp_add_sound_source!(
        simulation._cpp_simulation,
        source_args...,
        pq_type,
        _cmplx_split_and_store(fdpqv[])...
    )
end

import HDF5
import GLMakie
import LinearAlgebra

function order_points(
    cni::Vector{Int},
    cni_to_xyz::Vector{GLMakie.Point3f},
    dimension_to_eliminate::Int
)
    function segments_intersect(
        p1::Tuple{Float32,Float32},
        p2::Tuple{Float32,Float32},
        p3::Tuple{Float32,Float32},
        p4::Tuple{Float32,Float32}
    )
        function cross_prod(
            a::Tuple{Float32,Float32},
            b::Tuple{Float32,Float32},
            c::Tuple{Float32,Float32}
        )
            return (b[1] - a[1])*(c[2] - a[2]) - (b[2] - a[2])*(c[1] - a[1])
        end
        d1 = cross_prod(p3, p4, p1)
        d2 = cross_prod(p3, p4, p2)
        d3 = cross_prod(p1, p2, p3)
        d4 = cross_prod(p1, p2, p4)

        return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0))
    end

    dims = filter!(x -> x != dimension_to_eliminate, [1,2,3])
    a = (cni_to_xyz[cni[1]][dims[1]], cni_to_xyz[cni[1]][dims[2]])
    b = (cni_to_xyz[cni[2]][dims[1]], cni_to_xyz[cni[2]][dims[2]])
    c = (cni_to_xyz[cni[3]][dims[1]], cni_to_xyz[cni[3]][dims[2]])
    d = (cni_to_xyz[cni[4]][dims[1]], cni_to_xyz[cni[4]][dims[2]])

    if segments_intersect(a, c, b, d)
        return (cni[1], cni[2], cni[3], cni[4])
    end
    if segments_intersect(a, b, c, d)
        return (cni[1], cni[3], cni[2], cni[4])
    end
    if segments_intersect(a, d, b, c)
        return (cni[1], cni[2], cni[4], cni[3])
    end
    @assert false
end

function plot_pressure_field(h5_path::AbstractString)

    file = HDF5.h5open(h5_path, "r")

    ni_to_xyz = HDF5.read(file["/inputs/mesh/nodes"])
    ni_to_xyz = Float32.(ni_to_xyz)

    sei_to_ni = HDF5.read(file["/inputs/mesh/surface_elements"])
    sei_to_ni = convert(Matrix{Int}, sei_to_ni)
    sei_count::Int = size(sei_to_ni, 2)

    vei_to_ni = HDF5.read(file["/inputs/mesh/volume_elements"])
    vei_to_ni = convert(Matrix{Int}, vei_to_ni)
    vei_count::Int = size(vei_to_ni, 2)

    fi_to_freq = HDF5.read(file["/inputs/simulated_frequencies"])
    fi_to_freq = Float32.(fi_to_freq)
    fi_count::Int = length(fi_to_freq)

    # get segments
    segments::Dict{Tuple{Int,Int}, Vector{Int}} = Dict() # ni tuple to vei vec
    for vei in 1:vei_count
        for (eniv_1, eniv_2) in [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]
            ni_1 = vei_to_ni[eniv_1, vei]
            ni_2 = vei_to_ni[eniv_2, vei]
            if ni_1 > ni_2
                ni_1, ni_2 = ni_2, ni_1
            end
            if !haskey(segments, (ni_1, ni_2))
                segments[(ni_1, ni_2)] = [vei]
            else
                push!(segments[(ni_1, ni_2)], vei)
            end
        end
    end

    # compute the sliced surface for a single plane -- now returns fresh,
    # self-contained data instead of appending to shared arrays, so each
    # plane (x/y/z) can be computed, stored, and updated independently.
    function compute_plane_cut(dim::Int, position::Float32)
        cni_to_xyz_local::Vector{GLMakie.Point3f} = []
        cni_to_weights_local::Vector{Tuple{Int,Float32,Int,Float32}} = []
        faces_local::Vector{GLMakie.GLTriangleFace} = []
        cut_tets::Dict{Int, Vector{Int}} = Dict() # vei to cni vector

        for ((ni_1, ni_2), vei_vector) in segments
            dist_1::Float32 = ni_to_xyz[dim, ni_1] - position
            dist_2::Float32 = ni_to_xyz[dim, ni_2] - position
            if (dist_1 >= 0f0 && dist_2 < 0f0) || (dist_1 < 0f0 && dist_2 >= 0f0)
                t::Float32 = dist_1 / (dist_1 - dist_2)
                w_1::Float32 = 1 - t
                w_2::Float32 = t
                push!(cni_to_weights_local, (ni_1, w_1, ni_2, w_2))
                xyz_cut::GLMakie.Point3f = GLMakie.Point3f(
                    ni_to_xyz[:,ni_1] .+
                    t.*(ni_to_xyz[:,ni_2] .- ni_to_xyz[:,ni_1])
                )
                push!(cni_to_xyz_local, xyz_cut)
                cni = length(cni_to_xyz_local)
                for vei in vei_vector
                    if !haskey(cut_tets, vei)
                        cut_tets[vei] = [cni]
                    else
                        push!(cut_tets[vei], cni)
                    end
                end
            end
        end

        # loop over tets to create faces
        for (vei, cni_vector) in cut_tets
            if length(cni_vector) == 3
                push!(faces_local, GLMakie.GLTriangleFace(cni_vector...))
            elseif length(cni_vector) == 4
                function dist(cni1, cni2)
                    dims = filter!(x -> x != dim, [1,2,3])
                    return LinearAlgebra.norm(
                        cni_to_xyz_local[cni1][dims] .- cni_to_xyz_local[cni2][dims]
                    )
                end
                t = 1e-6
                (a, b, c, d) = cni_vector
                if (dist(a,b) < t)
                    push!(faces_local, GLMakie.GLTriangleFace(a, c, d))
                elseif (dist(a,c) < t) || (dist(b,c) < t)
                    push!(faces_local, GLMakie.GLTriangleFace(a, b, d))
                elseif (dist(a,d) < t) || (dist(b,d) < t) || (dist(c,d) < t)
                    push!(faces_local, GLMakie.GLTriangleFace(a, b, c))
                else
                    (cni_1, cni_2, cni_3, cni_4) =
                        order_points(cni_vector, cni_to_xyz_local, dim)
                    push!(faces_local, GLMakie.GLTriangleFace(cni_1, cni_3, cni_2))
                    push!(faces_local, GLMakie.GLTriangleFace(cni_1, cni_3, cni_4))
                end
            else
                @assert false
            end
        end

        return cni_to_xyz_local, faces_local, cni_to_weights_local
    end

    # color computation for a given plane's weights at a given frequency index
    function compute_colors(
        cni_to_weights_local::Vector{Tuple{Int,Float32,Int,Float32}},
        fi::Int
    )
        pressure = file["/results/pressure"][:,fi]
        colors_local = Vector{Float32}(undef, length(cni_to_weights_local))
        for cni in eachindex(cni_to_weights_local)
            (ni_1, w_1, ni_2, w_2) = cni_to_weights_local[cni]
            colors_local[cni] =
                abs(ComplexF32(pressure[ni_1])) * w_1 +
                abs(ComplexF32(pressure[ni_2])) * w_2
        end
        return colors_local
    end

    # get surface elements segments to plot boundaries
    sfc_segments_set::Set{Tuple{Int,Int}} = Set()
    for sei in 1:sei_count
        for (enis_1, enis_2) in [(1,2), (1,3), (2,3)]
            ni_1::Int = sei_to_ni[enis_1, sei]
            ni_2::Int = sei_to_ni[enis_2, sei]
            if ni_1 > ni_2
                ni_1, ni_2 = ni_2, ni_1
            end
            push!(sfc_segments_set, (ni_1, ni_2))
        end
    end
    sfc_segments::Vector{GLMakie.Point3f} = Vector{GLMakie.Point3f}(
        undef, 2*length(sfc_segments_set)
    )
    i::Int = 1
    for segment in sfc_segments_set
        ni1 = segment[1]
        ni2 = segment[2]
        xyz1 = GLMakie.Point3f(ni_to_xyz[:, ni1])
        xyz2 = GLMakie.Point3f(ni_to_xyz[:, ni2])
        sfc_segments[i] = xyz1
        i += 1
        sfc_segments[i] = xyz2
        i += 1
    end
    empty!(sfc_segments_set)

    # create figure
    fig = GLMakie.Figure()
    ax = GLMakie.Axis3(fig[1,1], aspect=:data, perspectiveness=0.5)
    ax.title = "Pressure field"

    # draw surface elements
    GLMakie.linesegments!(
        ax, sfc_segments;
        color = RGBA(0, 0, 0, 0.2),
        linewidth = 1,
        transparency = true,
    )

    GLMakie.Colorbar(
        fig[1,2];
        colormap = :rainbow1,
        colorrange = (0, 30),
        label = "Pressure",
    )

    # bounding box
    x_min = minimum(ni_to_xyz[1,:])
    x_max = maximum(ni_to_xyz[1,:])
    y_min = minimum(ni_to_xyz[2,:])
    y_max = maximum(ni_to_xyz[2,:])
    z_min = minimum(ni_to_xyz[3,:])
    z_max = maximum(ni_to_xyz[3,:])

    # per-plane mutable state
    fi::Ref{Int} = Ref(1)
    dim_names = ["x", "y", "z"]
    positions = [(x_min + x_max)/2, (y_min + y_max)/2, (z_min + z_max)/2]
    visible_flags = [true, true, true]
    mesh_objs::Vector{Any} = [nothing, nothing, nothing]
    weights_data::Vector{Vector{Tuple{Int,Float32,Int,Float32}}} = [[], [], []]

    function update_plane!(dim::Int)
        if mesh_objs[dim] !== nothing
            GLMakie.delete!(ax, mesh_objs[dim])
            mesh_objs[dim] = nothing
        end

        if !visible_flags[dim]
            weights_data[dim] = []
            return
        end

        cni_to_xyz_local, faces_local, cni_to_weights_local =
            compute_plane_cut(dim, positions[dim])
        weights_data[dim] = cni_to_weights_local

        if !isempty(cni_to_xyz_local)
            colors_local = compute_colors(cni_to_weights_local, fi[])
            mesh_objs[dim] = GLMakie.mesh!(
                ax, cni_to_xyz_local, faces_local;
                color = colors_local,
                colorrange = (0,30),
                shading = false,
                colormap = :rainbow1,
                transparency = false,
            )
        end
    end

    # frequency slider
    freq_slider_text = "Frequency: $(fi_to_freq[fi[]]) Hz (index $(fi[]))"
    slider = GLMakie.Slider(
        fig[2, 1],
        range = 1:fi_count,
        startvalue = fi[],
        horizontal = true,
    )
    label_slider = GLMakie.Label(
        fig[2, 2],
        freq_slider_text
    )
    on(slider.value) do fi_input
        fi[] = fi_input
        label_slider.text = freq_slider_text
        for dim in (1,2,3)
            if mesh_objs[dim] !== nothing
                mesh_objs[dim].color = compute_colors(weights_data[dim], fi[])
            end
        end
    end

    # --- Plane controls: label, position textbox, show/hide toggle ---
    limits = [(x_min, x_max), (y_min, y_max), (z_min, z_max)]
    controls_grid = fig[3, 1:2] = GLMakie.GridLayout()
    col = 1
    for dim in (1, 2, 3)
        GLMakie.Label(controls_grid[1, col], "$(dim_names[dim]):")
        col += 1

        # toggle
        toggle = GLMakie.Toggle(controls_grid[1, col], active = true)
        on(toggle.active) do active
            visible_flags[dim] = active
            update_plane!(dim)
        end
        col += 1

        # textbox
        (l, u) = limits[dim]
        textbox = GLMakie.Textbox(
            controls_grid[1, col],
            placeholder = first(string(l), 4) * " | " * first(string(u), 4),
            stored_string = first(string(positions[dim]), 5),
            validator = s -> begin
                v = tryparse(Float32, s)
                v !== nothing && l <= v <= u
            end,
            width = 100,
        )
        on(textbox.displayed_string) do s
            v = tryparse(Float32, s)
            if v !== nothing
                positions[dim] = v
                update_plane!(dim)
            end
        end
        col += 1
    end

    # initial draw for all three planes
    update_plane!(1)
    update_plane!(2)
    update_plane!(3)

    GLMakie.display(fig)

    closed = Ref(false)
    on(events(fig).window_open) do is_open
        if !is_open && !closed[]
            close(file)
            closed[] = true
        end
    end

    return fig
end
