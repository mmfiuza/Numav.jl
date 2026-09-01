# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export plot_pressure_field

import HDF5
import Statistics
import GLMakie
import LinearAlgebra

function _order_points(
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

function plot_pressure_field(
    file_path::AbstractString;
    db::Bool = false,
    colormap::Symbol = :rainbow1,
)
    file = HDF5.h5open(file_path, "r")

    get_sim_type(att::String) = HDF5.attrs(file["/simulation_type"])[att]
    if (
        get_sim_type("numerical_method") != "finite_element_method" ||
        get_sim_type("equation") != "helmholtz" ||
        get_sim_type("element_shape") != "tetrahedron" ||
        (   
            get_sim_type("element_order") != "linear" &&
            get_sim_type("element_order") != "quadratic"
        )
    )
        _throw_simulation_not_applicable()
    end

    ni_to_xyz = HDF5.read(file["/inputs/mesh/nodes"])
    ni_to_xyz = Float32.(ni_to_xyz)

    vei_to_ni = HDF5.read(file["/inputs/mesh/volume_elements"])
    vei_to_ni = convert(Matrix{Int}, vei_to_ni)
    vei_count::Int = size(vei_to_ni, 2)

    if get_sim_type("element_order") == "quadratic"
        @warn(
            "`plot_pressure_field` is not accurate for second order elements."*
            " The graph shown does not account for the non-vertex nodes"
        )
        vei_to_ni = vei_to_ni[1:4, :]
    end

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

    # bounding box
    (x_min::Float32, x_max::Float32) = extrema(ni_to_xyz[1,:])
    (y_min::Float32, y_max::Float32) = extrema(ni_to_xyz[2,:])
    (z_min::Float32, z_max::Float32) = extrema(ni_to_xyz[3,:])

    function pressure_calc_func(p::ComplexF32)
        if db
            return 20*log10( abs(p) / (sqrt(2)*20e-6) )
        else
            return abs(p)
        end
    end

    # pressure range
    p_all = pressure_calc_func.(ComplexF32.(file["/results/pressure"][:,:]))
    if db
        p_all = filter(isfinite, p_all) # remove all NaN
    end
    (p_min::Float32, p_max::Float32) = extrema(p_all)
    p_mean = Statistics.mean(p_all)
    p_std = Statistics.std(p_all)
    p_all = []

    # per-plane mutable state
    cni_to_xyz::Vector{Vector{GLMakie.Point3f}} = [[],[],[]]
    cni_to_weights::Vector{Vector{Tuple{Int,Float32,Int,Float32}}} = [[],[],[]]
    faces::Vector{Vector{GLMakie.GLTriangleFace}} = [[],[],[]]
    colors::Vector{Vector{Float32}} = [[],[],[]]
    fi::Int = 1
    dim_names::Vector{String} = ["yz", "xz", "xy"]
    positions::Vector{Float32} = [
        (x_min + x_max)/2, (y_min + y_max)/2, (z_min + z_max)/2
    ]
    color_range = (max(p_mean - 2.5*p_std, 0), p_mean + 2.5*p_std)
    visible_flags = [true, true, true]
    mesh_objs::Vector{Any} = [nothing, nothing, nothing]

    # compute the sliced surface for a single plane
    function compute_plane_cut(dim::Int, position::Float32)
        empty!(cni_to_xyz[dim])
        empty!(cni_to_weights[dim])
        empty!(faces[dim])

        cut_tets::Dict{Int, Vector{Int}} = Dict() # vei to cni vector
        for ((ni_1, ni_2), vei_vector) in segments
            dist1::Float32 = ni_to_xyz[dim, ni_1] - position
            dist2::Float32 = ni_to_xyz[dim, ni_2] - position
            if (dist1 >= 0f0 && dist2 < 0f0) || (dist1 < 0f0 && dist2 >= 0f0)
                t::Float32 = dist1 / (dist1 - dist2)
                w_1::Float32 = 1 - t
                w_2::Float32 = t
                push!(cni_to_weights[dim], (ni_1, w_1, ni_2, w_2))
                xyz_cut::GLMakie.Point3f = GLMakie.Point3f(
                    ni_to_xyz[:,ni_1] .+
                    t.*(ni_to_xyz[:,ni_2] .- ni_to_xyz[:,ni_1])
                )
                push!(cni_to_xyz[dim], xyz_cut)
                cni = length(cni_to_xyz[dim])
                for vei in vei_vector
                    if !haskey(cut_tets, vei)
                        cut_tets[vei] = [cni]
                    else
                        push!(cut_tets[vei], cni)
                    end
                end
            end
        end
        empty!(colors[dim])
        colors[dim] = Vector{Float32}(undef, length(cni_to_xyz[dim]))

        # loop over tets to create faces
        for (vei, cni_vector) in cut_tets
            if length(cni_vector) == 3
                push!(faces[dim], GLMakie.GLTriangleFace(cni_vector...))
            elseif length(cni_vector) == 4
                function dist(cni1, cni2)
                    dims = filter!(x -> x != dim, [1,2,3])
                    return LinearAlgebra.norm(
                        cni_to_xyz[dim][cni1][dims] .-
                        cni_to_xyz[dim][cni2][dims]
                    )
                end
                t = 1e-6
                (a, b, c, d) = cni_vector
                if (dist(a,b) < t)
                    push!(faces[dim], GLMakie.GLTriangleFace(a, c, d))
                elseif (dist(a,c) < t) || (dist(b,c) < t)
                    push!(faces[dim], GLMakie.GLTriangleFace(a, b, d))
                elseif (dist(a,d) < t) || (dist(b,d) < t) || (dist(c,d) < t)
                    push!(faces[dim], GLMakie.GLTriangleFace(a, b, c))
                else
                    (cni_1, cni_2, cni_3, cni_4) =
                        _order_points(cni_vector, cni_to_xyz[dim], dim)
                    push!(faces[dim], GLMakie.GLTriangleFace(cni_1, cni_3, cni_2))
                    push!(faces[dim], GLMakie.GLTriangleFace(cni_1, cni_3, cni_4))
                end
            else
                @assert false
            end
        end
    end

    # color computation for a given plane's weights at a given frequency index
    function compute_colors(dim::Int, fi::Int)
        pressure = file["/results/pressure"][:,fi]
        for cni in eachindex(cni_to_weights[dim])
            (ni_1, w_1, ni_2, w_2) = cni_to_weights[dim][cni]
            colors[dim][cni] =
                pressure_calc_func(
                    ComplexF32(pressure[ni_1]) * w_1 +
                    ComplexF32(pressure[ni_2]) * w_2
                )
        end
        return colors[dim]
    end

    # find the boundary faces of the volume mesh
    face_count::Dict{NTuple{3,Int}, Int} = Dict()
    for vei in 1:vei_count
        for (eniv_1, eniv_2, eniv_3) in [(1,2,3), (1,2,4), (1,3,4), (2,3,4)]
            ni_1 = vei_to_ni[eniv_1, vei]
            ni_2 = vei_to_ni[eniv_2, vei]
            ni_3 = vei_to_ni[eniv_3, vei]
            face_key = Tuple(sort([ni_1, ni_2, ni_3]))
            face_count[face_key] = get(face_count, face_key, 0) + 1
        end
    end

    # get surface elements segments to plot boundaries
    sfc_segments_set::Set{Tuple{Int,Int}} = Set()
    for (face_key, count) in face_count
        if count == 1
            (ni_1, ni_2, ni_3) = face_key
            for (eni_1, eni_2) in [(ni_1,ni_2), (ni_1,ni_3), (ni_2,ni_3)]
                if eni_1 > eni_2
                    eni_1, eni_2 = eni_2, eni_1
                end
                push!(sfc_segments_set, (eni_1, eni_2))
            end
        end
    end
    empty!(face_count)
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
    GLMakie.activate!(title = "Pressure field")
    fig = GLMakie.Figure()
    ax = GLMakie.Axis3(fig[1,1], aspect=:data, perspectiveness=0.5)
    ax.title = ""

    # draw surface elements
    GLMakie.linesegments!(
        ax, sfc_segments;
        color = GLMakie.RGBA(0, 0, 0, 0.2),
        linewidth = 1,
        transparency = true,
    )

    color_bar_label = if db
        "Sound pressure level (dB, ref. 20 \u03bcPa)"
    else
        "Sound pressure amplitude (Pa)"
    end

    # color bar
    colorbar = GLMakie.Colorbar(
        fig[1,2];
        colormap = colormap,
        colorrange = color_range,
        label = color_bar_label,
        ticks = GLMakie.LinearTicks(10),
    )
    range_slider = GLMakie.IntervalSlider(
        fig[1,3],
        range = LinRange(max(p_min, 0), (p_mean + 5*p_std), 256),
        startvalues = color_range,
        horizontal = false,
    )
    GLMakie.on(range_slider.interval) do interval
        color_range = interval
        colorbar.colorrange = color_range
        for dim in (1,2,3)
            if mesh_objs[dim] !== nothing
                mesh_objs[dim].colorrange = color_range
            end
        end
    end

    function update_plane(dim::Int)
        if mesh_objs[dim] !== nothing
            GLMakie.delete!(ax, mesh_objs[dim])
            mesh_objs[dim] = nothing
        end
        if !visible_flags[dim]
            return
        end

        compute_plane_cut(dim, positions[dim])
        if !isempty(cni_to_xyz[dim])
            compute_colors(dim, fi)
            mesh_objs[dim] = GLMakie.mesh!(
                ax, cni_to_xyz[dim], faces[dim];
                color = colors[dim],
                colorrange = color_range,
                shading = false,
                colormap = colormap,
                transparency = false,
            )
        end
    end

    # frequency slider
    freq_slider_text(fi::Int) = "Frequency: $(fi_to_freq[fi]) Hz (index $(fi))"
    slider = GLMakie.Slider(
        fig[2, 1],
        range = 1:fi_count,
        startvalue = fi,
        horizontal = true,
    )
    label_slider = GLMakie.Label(
        fig[2, 2],
        freq_slider_text(1)
    )
    GLMakie.on(slider.value) do fi_input
        fi = fi_input
        label_slider.text = freq_slider_text(fi)
        for dim in (1,2,3)
            if mesh_objs[dim] !== nothing
                compute_colors(dim, fi)
                mesh_objs[dim].color = colors[dim]
            end
        end
    end

    # plane controls: label, position textbox, show/hide toggle
    limits = [(x_min, x_max), (y_min, y_max), (z_min, z_max)]
    controls_grid = fig[3, 1:2] = GLMakie.GridLayout()
    col = 1
    for dim in (3, 2, 1)
        GLMakie.Label(controls_grid[1, col], "$(dim_names[dim]):")
        col += 1

        # toggle
        toggle = GLMakie.Toggle(controls_grid[1, col], active = true)
        GLMakie.on(toggle.active) do active
            visible_flags[dim] = active
            update_plane(dim)
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
        GLMakie.on(textbox.displayed_string) do s
            v = tryparse(Float32, s)
            if v !== nothing
                positions[dim] = v
                update_plane(dim)
            end
        end
        col += 1
    end

    # initial draw for all three planes
    update_plane(1)
    update_plane(2)
    update_plane(3)

    GLMakie.display(fig)

    closed = Ref(false)
    GLMakie.on(GLMakie.events(fig).window_open) do is_open
        if !is_open && !closed[]
            close(file)
            closed[] = true
        end
    end

    return fig
end
