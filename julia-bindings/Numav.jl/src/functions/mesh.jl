# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file has functions applicable to Simulation subtypes that have a mesh.

export load_mesh!

const SimulationWithMesh = Union{SimulationFemHelmholtz}

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
        # TODO: change 1 to undef
        s._sei_to_ni = ones(_enis_count(s), s._sei_count)
        s._vei_to_ni = ones(_eniv_count(s), s._vei_count)
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
                s._sei_to_espg[sei] = espg
                s._sei_to_ni[1:3, sei] = [
                    parse(Int, strip(line[25:32])),
                    parse(Int, strip(line[33:40])),
                    parse(Int, strip(line[41:48])),
                ]
                sei += 1

            elseif startswith(line, "CTETRA  ")
                evpg = parse(Int, strip(line[17:24]))
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

function load_mesh!(
    s::SimulationWithMesh,
    path_to_mesh::AbstractString
)   
    (_, ext) = splitext(path_to_mesh)
    if ext == ".bdf" || ext == ".nas"
        _load_bdf!(s, path_to_mesh)
    else
        if ext == ""
            throw(ArgumentError("No file format on file name", ))
        else
            throw(ArgumentError("Unrecognized file format: `$ext`"))
        end
    end
    _cpp_load_mesh(
        s._cpp_simulation,
        vec(s._ni_to_xyz),
        UInt64.(vec(s._sei_to_ni) .- 1),
        UInt64.(vec(s._vei_to_ni) .- 1),
        UInt64.(s._sei_to_espg),
        UInt64.(s._vei_to_evpg),
        UInt64(s._ni_count),
        UInt64(s._sei_count),
        UInt64(s._vei_count)
    )
end
