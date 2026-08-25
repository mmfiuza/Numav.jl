# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file has functions applicable to Simulation subtypes that have a mesh.

export load_mesh!

const SimulationWithMesh = Union{SimulationFemHelmholtz}

function _check_if_mesh_is_defined(s::SimulationWithMesh)
    if !s._is_mesh_defined
        error("Mesh not defined. Call load_mesh! to do so.")
    end
end

function _load_bdf!(
    s::SimulationWithMesh,
    path_to_mesh::AbstractString,
)
    if !isfile(path_to_mesh)
        error("Could not open file: $path_to_mesh")
    end

    open(path_to_mesh, "r") do file
        # first pass: count lines by type
        s._ni_count = 0
        s._sei_count = 0
        s._vei_count = 0
        for line in eachline(file)
            if     startswith(line, "GRID    ")
                s._ni_count += 1
            elseif startswith(line, "CTRIA3  ")
                s._sei_count += 1
            elseif startswith(line, "CTETRA  ")
                s._vei_count += 1
            elseif startswith(line, "ENDDATA ")
                break
            end
        end

        s._ni_to_xyz = Matrix{Float64}(undef, 3, s._ni_count)
        s._sei_to_ni = Matrix{Int}(undef, _enis_count(s), s._sei_count)
        s._vei_to_ni = Matrix{Int}(undef, _eniv_count(s), s._vei_count)
        s._sei_to_espg = Vector{Int}(undef, s._sei_count)
        s._vei_to_evpg = Vector{Int}(undef, s._vei_count)

        sei::Int = 1
        vei::Int = 1

        # second pass: parse data
        seekstart(file)
        for line in eachline(file)
            if startswith(line, "GRID    ")
                ni = Int(parse(Float64, strip(line[9:16])))
                s._ni_to_xyz[:, ni] = [
                    parse(Float64, strip(line[25:32])),
                    parse(Float64, strip(line[33:40])),
                    parse(Float64, strip(line[41:48])),
                ]

            elseif startswith(line, "CTRIA3  ")
                espg = parse(Int, strip(line[17:24]))
                push!(s._existing_espg, espg)
                s._sei_to_espg[sei] = espg
                s._sei_to_ni[1:3, sei] = [
                    parse(Int, strip(line[25:32])),
                    parse(Int, strip(line[33:40])),
                    parse(Int, strip(line[41:48])),
                ]
                sei += 1

            elseif startswith(line, "CTETRA  ")
                evpg = parse(Int, strip(line[17:24]))
                push!(s._existing_evpg, evpg)
                s._vei_to_evpg[vei] = evpg
                s._vei_to_ni[1:4, vei] = [
                    parse(Int, strip(line[25:32])),
                    parse(Int, strip(line[33:40])),
                    parse(Int, strip(line[41:48])),
                    parse(Int, strip(line[49:56])),
                ]
                vei += 1

            elseif startswith(line, "ENDDATA ")
                break
            end
        end
    end
end

function _generate_non_vtx_nodes!(s::SimulationFemHelmholtz)
    # 1-based vertex-pair tables (translated from the 0-based C++ arrays)
    VTX_PAIRS_VOL = ((1,2), (1,3), (1,4), (2,3), (2,4), (3,4))
    VTX_PAIRS_SFC = ((1,2), (1,3), (2,3))

    ENIS_COUNT_LIN = _enis_count(ElementOrder_LINEAR)
    ENIV_COUNT_LIN = _eniv_count(ElementOrder_LINEAR)

    idxs_extra_nodes = Dict{Tuple{Int,Int}, Int}()

    # first pass: count extra nodes and save idx tuples
    is_extra_node = falses(length(VTX_PAIRS_VOL), s._vei_count)
    for vei in 1:s._vei_count
        for i in eachindex(VTX_PAIRS_VOL)
            a, b = VTX_PAIRS_VOL[i]
            n1 = s._vei_to_ni[a, vei]
            n2 = s._vei_to_ni[b, vei]
            tup = n1 < n2 ? (n1, n2) : (n2, n1)

            if !haskey(idxs_extra_nodes, tup)
                is_extra_node[i, vei] = true
                s._ni_count += 1
                s._vei_to_ni[ENIV_COUNT_LIN + i, vei] = s._ni_count
                idxs_extra_nodes[tup] = s._ni_count
            else
                is_extra_node[i, vei] = false
                s._vei_to_ni[ENIV_COUNT_LIN + i, vei] = idxs_extra_nodes[tup]
            end
        end
    end

    # grow _ni_to_xyz to fit the newly counted nodes
    ni_to_xyz_new = zeros(Float64, 3, s._ni_count)
    ni_to_xyz_new[:, 1:size(s._ni_to_xyz, 2)] .= s._ni_to_xyz
    s._ni_to_xyz = ni_to_xyz_new

    # second pass: create the extra nodes and assign to volume elements
    for vei in 1:s._vei_count
        for i in eachindex(VTX_PAIRS_VOL)
            is_extra_node[i, vei] || continue

            a, b = VTX_PAIRS_VOL[i]
            n1 = s._vei_to_ni[a, vei]
            n2 = s._vei_to_ni[b, vei]
            tup = n1 < n2 ? (n1, n2) : (n2, n1)

            x = (s._ni_to_xyz[1, tup[1]] + s._ni_to_xyz[1, tup[2]]) / 2
            y = (s._ni_to_xyz[2, tup[1]] + s._ni_to_xyz[2, tup[2]]) / 2
            z = (s._ni_to_xyz[3, tup[1]] + s._ni_to_xyz[3, tup[2]]) / 2

            idx_extra_node = idxs_extra_nodes[tup]
            s._ni_to_xyz[:, idx_extra_node] = [x, y, z]
        end
    end

    # third pass: assign nodes to surface elements
    for sei in 1:s._sei_count
        for i in eachindex(VTX_PAIRS_SFC)
            a, b = VTX_PAIRS_SFC[i]
            n1 = s._sei_to_ni[a, sei]
            n2 = s._sei_to_ni[b, sei]
            tup = n1 < n2 ? (n1, n2) : (n2, n1)
            s._sei_to_ni[ENIS_COUNT_LIN + i, sei] = idxs_extra_nodes[tup]
        end
    end
end

function load_mesh!(
    s::SimulationWithMesh,
    path_to_mesh::AbstractString
)   
    if s._is_mesh_defined
        error("Mesh already defined")
    end
    (_, ext) = splitext(path_to_mesh)
    if ext == ".bdf" || ext == ".nas"
        _load_bdf!(s, path_to_mesh)
    else
        if ext == ""
            throw(ArgumentError("No file format on file name"))
        else
            throw(ArgumentError("Unrecognized file format: `$ext`"))
        end
    end
    if (_is_quadratic(s))
        _generate_non_vtx_nodes!(s)
    end
    s._is_mesh_defined = true
end
