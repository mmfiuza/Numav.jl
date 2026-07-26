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

# # --- Quadratic root helper (for second‑order edges) ---
# function find_roots_quadratic(a, b, c, tol)
#     roots = Float64[]
#     if abs(a) < 1e-12 && abs(b) > 1e-12
#         t = -c / b
#         if -tol <= t <= 1 + tol
#             t_clamped = clamp(t, 0.0, 1.0)
#             # check if the function value is really zero
#             if abs(a*t_clamped^2 + b*t_clamped + c) < 100*tol
#                 push!(roots, t_clamped)
#             end
#         end
#     else
#         D = b^2 - 4*a*c
#         if D >= -tol
#             D = max(D, 0.0)
#             sqrtD = sqrt(D)
#             t1 = (-b + sqrtD) / (2a)
#             t2 = (-b - sqrtD) / (2a)
#             for t in (t1, t2)
#                 if -tol <= t <= 1 + tol
#                     t_clamped = clamp(t, 0.0, 1.0)
#                     if abs(a*t_clamped^2 + b*t_clamped + c) < 100*tol
#                         push!(roots, t_clamped)
#                     end
#                 end
#             end
#         end
#     end
#     return unique(roots)  # remove duplicates caused by clamping
# end

function plot_pressure_field(
    h5_path::AbstractString;
    plane::Symbol = :xy,
    position::Real = 0.0,
    freq::Union{Real, Nothing} = nothing
)
    tolerance::Float32 = 1e-6

    HDF5.h5open(h5_path, "r") do file
        # --- Read mesh and data ---
        ni_to_xyz = HDF5.read(file["/inputs/mesh/nodes"])
        ni_to_xyz = Float32.(ni_to_xyz)

        sei_to_ni = HDF5.read(file["/inputs/mesh/surface_elements"])
        sei_to_ni = convert(Matrix{Int}, sei_to_ni)
        sei_count = size(sei_to_ni, 2)

        vei_to_ni = HDF5.read(file["/inputs/mesh/volume_elements"])
        vei_to_ni = convert(Matrix{Int}, vei_to_ni)
        vei_count = size(vei_to_ni, 2)

        pressure = HDF5.read(file["/results/pressure"])
        pressure = Float32.(abs.(pressure))

        fi_to_freq = HDF5.read(file["/inputs/simulated_frequencies"])
        fi_to_freq = Float32.(fi_to_freq)
        fi_count = length(fi_to_freq)

        # Determine initial frequency
        fi_initial = 1
        freq_initial = fi_to_freq[fi_initial]

        # Plane definition
        dim = if plane == :xy
            3
        elseif plane == :yz
            1
        elseif plane == :xz
            2
        else
            error("plane must be :xy, :yz or :xz")
        end
        dist(ni) = ni_to_xyz[dim,ni] - position

        # --- Detect element order ---
        is_order_linear = size(vei_to_ni, 1) == 4

        # --- Slice volume mesh ---
        verts_all = GLMakie.Point3f[]
        faces_all = GLMakie.GLTriangleFace[]
        vertex_mappings = Vector{Vector{Tuple{Int,Float32}}}()

        for vei in 1:vei_count
            eniv_to_ni = @view vei_to_ni[:, vei]   # all node indices (4 or 10)

            plane_points = Tuple{Float32,Float32,Float32}[]
            plane_contribs = Vector{Tuple{Int,Float32}}[]

            if is_order_linear
                eniv_to_dist = [dist(eniv_to_ni[eniv]) for eniv in 1:4]
                pos_mask = eniv_to_dist .> tolerance
                neg_mask = eniv_to_dist .< -tolerance
                zero_mask = abs.(eniv_to_dist) .<= tolerance

                if !((any(pos_mask) && any(neg_mask)) || (count(zero_mask) > 2))
                    continue
                end

                tet_coords = ni_to_xyz[:, eniv_to_ni]
                eei_to_eniv = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]
                for (eniv_1, eniv_2) in eei_to_eniv
                    dist_1 = eniv_to_dist[eniv_1]
                    dist_2 = eniv_to_dist[eniv_2]
                    if abs(dist_1) <= tolerance && abs(dist_2) <= tolerance
                        push!(plane_points, Tuple(tet_coords[:,eniv_1]), Tuple(tet_coords[:,eniv_2]))
                        push!(plane_contribs, [(eniv_to_ni[eniv_1], 1.0f0)], [(eniv_to_ni[eniv_2], 1.0f0)])
                    elseif abs(dist_1) <= tolerance
                        push!(plane_points, Tuple(tet_coords[:,eniv_1]))
                        push!(plane_contribs, [(eniv_to_ni[eniv_1], 1.0f0)])
                    elseif abs(dist_2) <= tolerance
                        push!(plane_points, Tuple(tet_coords[:,eniv_2]))
                        push!(plane_contribs, [(eniv_to_ni[eniv_2], 1.0f0)])
                    elseif (dist_1 > 0) != (dist_2 > 0)
                        t = dist_1 / (dist_1 - dist_2)
                        xyz_intersect = tet_coords[:,eniv_1] .+ t .* (tet_coords[:,eniv_2] .- tet_coords[:,eniv_1])
                        push!(plane_points, Tuple(xyz_intersect))
                        push!(plane_contribs, [(eniv_to_ni[eniv_1], 1.0f0 - t), (eniv_to_ni[eniv_2], t)])
                    end
                end
            else
                # # Quadratic tetrahedron slicing
                # eei_to_eniv = [
                #     (1,2,5),
                #     (1,3,6),
                #     (1,4,7),
                #     (2,3,8),
                #     (2,4,9),
                #     (3,4,10)
                # ]
                # for (eniv_1, eniv_2, eniv_m) in eei_to_eniv
                #     ni_1 = eniv_to_ni[eniv_1]
                #     ni_2 = eniv_to_ni[eniv_2]
                #     ni_m = eniv_to_ni[eniv_m]

                #     # slicing coordinate for each node
                #     p_1 = ni_to_xyz[dim, ni_1]
                #     p_2 = ni_to_xyz[dim, ni_2]
                #     p_m = ni_to_xyz[dim, ni_m]

                #     # coefficients of quadratic d(t) = a t^2 + b t + (ci - position)
                #     a = 2p_1 + 2p_2 - 4p_m
                #     b = -3p_1 - p_2 + 4p_m
                #     c = p_1 - position

                #     roots = find_roots_quadratic(a, b, c, tolerance)
                #     for t in roots
                #         # quadratic shape functions along the edge
                #         w_1 = (1 - t) * (1 - 2t)
                #         w_2 = t * (2t - 1)
                #         w_m = 4 * t * (1 - t)

                #         # intersection point
                #         x =
                #             ni_to_xyz[1,ni_1]*w_1 +
                #             ni_to_xyz[1,ni_2]*w_2 +
                #             ni_to_xyz[1,ni_m]*w_m
                #         y =
                #             ni_to_xyz[2,ni_1]*w_1 +
                #             ni_to_xyz[2,ni_2]*w_2 +
                #             ni_to_xyz[2,ni_m]*w_m
                #         z =
                #             ni_to_xyz[3,ni_1]*w_1 +
                #             ni_to_xyz[3,ni_2]*w_2 +
                #             ni_to_xyz[3,ni_m]*w_m
                #         push!(plane_points, (x, y, z))

                #         # pressure mapping (only non‑negligible weights)
                #         contribs = Tuple{Int,Float32}[]
                #         if abs(w_1) > tolerance
                #             push!(contribs, (ni_1, w_1))
                #         end
                #         if abs(w_2) > tolerance
                #             push!(contribs, (ni_2, w_2))
                #         end
                #         if abs(w_m) > tolerance
                #             push!(contribs, (ni_m, w_m))
                #         end
                #         push!(plane_contribs, contribs)
                #     end
                # end
            end

            # --- Remove duplicate intersection points ---
            unique_pts = Tuple{Float32,Float32,Float32}[]
            unique_maps = Vector{Tuple{Int,Float32}}[]
            for (k, pt) in enumerate(plane_points)
                is_new = true
                for existing in unique_pts
                    if LinearAlgebra.norm(collect(pt) .- collect(existing)) < 10*tolerance
                        is_new = false
                        break
                    end
                end
                if is_new
                    push!(unique_pts, pt)
                    push!(unique_maps, plane_contribs[k])
                end
            end

            n_pts = length(unique_pts)
            if n_pts < 3
                continue
            end

            # Order points around centroid
            centroid = sum(collect.(unique_pts)) ./ n_pts
            if plane == :xy
                u = Float32.([1, 0, 0])
                v = Float32.([0, 1, 0])
            elseif plane == :yz
                u = Float32.([0, 1, 0])
                v = Float32.([0, 0, 1])
            else
                u = Float32.([1, 0, 0])
                v = Float32.([0, 0, 1])
            end
            rel_pos = [collect(pt) .- centroid for pt in unique_pts]
            angles = [atan(LinearAlgebra.dot(r, v), LinearAlgebra.dot(r, u)) for r in rel_pos]
            perm = sortperm(angles)
            ordered_pts = unique_pts[perm]
            ordered_maps = unique_maps[perm]

            v_start = length(verts_all) + 1
            for (pt, mp) in zip(ordered_pts, ordered_maps)
                push!(verts_all, GLMakie.Point3f(pt))
                push!(vertex_mappings, mp)
            end
            for k in 2:(n_pts-1)
                push!(faces_all, GLMakie.GLTriangleFace(v_start, v_start+k-1, v_start+k))
            end
        end

        # color computation
        n = length(vertex_mappings)
        colors = Vector{Float32}(undef, n)
        function compute_colors(fi::Int)
            @inbounds for i in 1:n
                colors[i] = 0.0f0 + 0.0f0*im
                for (ni, w) in vertex_mappings[i]
                    colors[i] += w * pressure[ni,fi]
                end
            end
            return colors
        end

        # create figure
        fig = GLMakie.Figure()
        ax = GLMakie.Axis3(fig[1,1], aspect=:data, perspectiveness=0.5)
        ax.title = "Pressure field"

        # draw surface elements
        segments = GLMakie.Point3f[]
        for sei in 1:sei_count
            ni1 = sei_to_ni[1, sei]
            ni2 = sei_to_ni[2, sei]
            ni3 = sei_to_ni[3, sei]
            xyz1 = GLMakie.Point3f(ni_to_xyz[:, ni1])
            xyz2 = GLMakie.Point3f(ni_to_xyz[:, ni2])
            xyz3 = GLMakie.Point3f(ni_to_xyz[:, ni3])
            push!(segments, xyz1, xyz2, xyz2, xyz3, xyz3, xyz1)
        end
        GLMakie.linesegments!(
            ax, segments;
            color = RGBA(0, 0, 0, 0.2),
            linewidth = 1,
            transparency = true,
        )

        mesh_obj = nothing
        if !isempty(verts_all)
            mesh_obj = GLMakie.mesh!(
                ax, verts_all, faces_all;
                color = compute_colors(fi_initial),
                # colorrange = global_limits,
                shading = false,
                colormap = :rainbow1,
                transparency = false,
            )
        end

        if mesh_obj !== nothing
            GLMakie.Colorbar(fig[1,2], mesh_obj, label="Pressure")
        end

        # --- Slider ---
        slider = GLMakie.Slider(
            fig[2, 1],
            range = 1:fi_count,
            startvalue = fi_initial,
            horizontal = true,
        )
        label_slider = GLMakie.Label(
            fig[2, 2],
            "Frequency: $(freq_initial) Hz (index $fi_initial)"
        )

        on(slider.value) do fi
            freq = fi_to_freq[fi]
            label_slider.text = "Frequency: $(freq) Hz (index $fi)"
            if mesh_obj !== nothing
                mesh_obj.color = compute_colors(fi)
            end
        end

        GLMakie.display(fig)
        return fig
    end
end
