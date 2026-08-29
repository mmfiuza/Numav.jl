# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file has functions applicable for any Simulation subtype

export
    create_simulation,
    set_result_export_path!,
    run!

function create_simulation(;
    numerical_method::Option,
    equation::Option,
    element_shape::Option,
    element_order::Option
)
    if numerical_method === Fem && equation === Helmholtz

        if element_shape === Tetrahedron
            element_shape_type = ElementShape_TETRAHEDRON
        else
            throw(ArgumentError("Invalid `element_shape` option"))
        end

        if element_order === Linear
            element_order_type = ElementOrder_LINEAR
        elseif element_order === Quadratic
            element_order_type = ElementOrder_QUADRATIC
        else
            throw(ArgumentError("Invalid `element_order` option"))
        end

        return SimulationFemHelmholtz{
            element_shape_type,
            element_order_type,
            _enis_count(element_order_type),
            _eniv_count(element_order_type),
        }()
    end
    throw(ArgumentError("Invalid options"))
end

function set_result_export_path!(
    s::Simulation,
    path_to_hdf5_file::AbstractString
)
    _check_if_did_run(s)
    if !isempty(s._hdf5_file_path)
        error("Result export path is already defined.")
    end
    s._hdf5_file_path = path_to_hdf5_file
end

function _organize_volume_physical_group_data!(s::SimulationFemHelmholtz)
    s._vei_to_ivpg = Vector{Int}(undef, s._vei_count)
    for vei in 1:s._vei_count
        s._vei_to_ivpg[vei] = s._evpg_ivpg_bimap[s._vei_to_evpg[vei]]
    end
end

function _organize_impedance_physical_group_data!(s::SimulationFemHelmholtz)
    s._isei_count = 0
    for sei in 1:s._sei_count
        if haskey(s._espg_ispgi_bimap, s._sei_to_espg[sei])
            s._isei_count += 1
        end
    end

    s._isei_to_ni = Matrix{Int}(undef, _enis_count(s), s._isei_count)
    s._isei_to_ispgi = Vector{Int}(undef, s._isei_count)
    isei = 0
    for sei in 1:s._sei_count
        if haskey(s._espg_ispgi_bimap, s._sei_to_espg[sei])
            isei += 1
            s._isei_to_ni[:, isei] = s._sei_to_ni[:, sei]
            s._isei_to_ispgi[isei] = s._espg_ispgi_bimap[s._sei_to_espg[sei]]
        end
    end
end

function _organize_velocity_physical_group_data!(s::SimulationFemHelmholtz)
    s._vsei_count = 0
    for sei in 1:s._sei_count
        if haskey(s._espg_ispgv_bimap, s._sei_to_espg[sei])
            s._vsei_count += 1
        end
    end

    s._vsei_to_ni = Matrix{Int}(undef, _enis_count(s), s._vsei_count)
    s._vsei_to_ispgv = Vector{Int}(undef, s._vsei_count)
    vsei = 0
    for sei in 1:s._sei_count
        if haskey(s._espg_ispgv_bimap, s._sei_to_espg[sei])
            vsei += 1
            s._vsei_to_ni[:, vsei] = s._sei_to_ni[:, sei]
            s._vsei_to_ispgv[vsei] = s._espg_ispgv_bimap[s._sei_to_espg[sei]]
        end
    end
end

function _organize_pressure_physical_group_data!(s::SimulationFemHelmholtz)
    s._psei_count = 0
    for sei in 1:s._sei_count
        if haskey(s._espg_ispgp_bimap, s._sei_to_espg[sei])
            s._psei_count += 1
        end
    end

    s._psei_to_ni = Matrix{Int}(undef, _enis_count(s), s._psei_count)
    s._psei_to_ispgp = Vector{Int}(undef, s._psei_count)
    psei = 0
    for sei in 1:s._sei_count
        if haskey(s._espg_ispgp_bimap, s._sei_to_espg[sei])
            psei += 1
            s._psei_to_ni[:, psei] = s._sei_to_ni[:, sei]
            s._psei_to_ispgp[psei] = s._espg_ispgp_bimap[s._sei_to_espg[sei]]
        end
    end

    # create vector of ni sets
    ni_sets::Vector{Set{Int}} = [Set{Int}() for _ in 1:s._ispgp_count]
    for psei in 1:s._psei_count
        ni = s._psei_to_ni[:, psei]
        ispgp::Int = s._psei_to_ispgp[psei]
        push!(ni_sets[ispgp], ni...)
    end
    for ppi in 1:s._ppi_count
        ni::Int = s._ppi_to_ni[ppi]
        push!(ni_sets, Set{Int}(ni))
    end

    ni_to_set_idxs = Dict{Int,Vector{Int}}()
    for (idx, set) in enumerate(ni_sets)
        for ni in set
            push!(get!(ni_to_set_idxs, ni, Int[]), idx)
        end
    end

    # group nis by their exact "membership signature" (which sets contain them)
    signature_to_ni = Dict{Vector{Int},Vector{Int}}()
    for (ni, idxs) in ni_to_set_idxs
        sig = sort(idxs)
        push!(get!(signature_to_ni, sig, Int[]), ni)
    end
    sorted_sigs = sort(collect(keys(signature_to_ni)))

    s._pvi_count = length(sorted_sigs)
    s._pvi_to_pressure_func = Vector{Function}(undef, s._pvi_count)
    s._pvi_to_pni_count = Vector{Int}(undef, s._pvi_count)
    s._pni_to_ni = Int[]

    pvi = 0
    for sig in sorted_sigs
        pvi += 1
        nis = signature_to_ni[sig]

        funcs::Vector{Function} = Function[]
        for idx in sig
            if idx <= s._ispgp_count
                push!(funcs, s._ispgp_to_pressure_func[idx])
            else
                ppi = idx - s._ispgp_count
                push!(funcs, s._ppi_to_pressure_func[ppi])
            end
        end

        n_funcs = length(funcs)
        s._pvi_to_pressure_func[pvi] = (freq) ->
            sum(f(freq) for f in funcs) / n_funcs

        s._pvi_to_pni_count[pvi] = length(nis)

        append!(s._pni_to_ni, sort(nis))
    end
    s._pni_count = length(s._pni_to_ni)

    return nothing
end

function _write_pq_matrices(s::SimulationFemHelmholtz)
    # _ivpg_to_density_values
    s._ivpg_to_density_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._ivpg_count, s._fi_count)
    for fi in 1:s._fi_count
        for ivpg in 1:s._ivpg_count
            s._ivpg_to_density_values[ivpg, fi] = 
                s._ivpg_to_density_func[ivpg](s._fi_to_freq[fi])
        end
    end
    
    # _ivpg_to_soundspeed_values
    s._ivpg_to_soundspeed_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._ivpg_count, s._fi_count)
    for ivpg in 1:s._ivpg_count
        for fi in 1:s._fi_count
            s._ivpg_to_soundspeed_values[ivpg, fi] =
                s._ivpg_to_soundspeed_func[ivpg](s._fi_to_freq[fi])
        end
    end

    # _ispgi_to_impedance_values
    s._ispgi_to_impedance_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._ispgi_count, s._fi_count)
    for ispgi in 1:s._ispgi_count
        for fi in 1:s._fi_count
            s._ispgi_to_impedance_values[ispgi, fi] =
                s._ispgi_to_impedance_func[ispgi](s._fi_to_freq[fi])
        end
    end

    # _ispgv_to_velocity_values
    s._ispgv_to_velocity_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._ispgv_count, s._fi_count)
    for ispgv in 1:s._ispgv_count
        for fi in 1:s._fi_count
            s._ispgv_to_velocity_values[ispgv, fi] =
                s._ispgv_to_velocity_func[ispgv](s._fi_to_freq[fi])
        end
    end

    # _vpi_to_volvel_values
    s._vpi_to_volvel_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._vpi_count, s._fi_count)
    for vpi in 1:s._vpi_count
        for fi in 1:s._fi_count
            s._vpi_to_volvel_values[vpi, fi] =
                s._vpi_to_volvel_func[vpi](s._fi_to_freq[fi])
        end
    end

    # _pvi_to_pressure_values
    s._pvi_to_pressure_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._pvi_count, s._fi_count)
    for pvi in 1:s._pvi_count
        for fi in 1:s._fi_count
            s._pvi_to_pressure_values[pvi, fi] =
                s._pvi_to_pressure_func[pvi](s._fi_to_freq[fi])
        end
    end
end


function run!(s::Simulation)
    _check_if_it_can_run(s)

    _organize_volume_physical_group_data!(s)
    _organize_impedance_physical_group_data!(s)
    _organize_velocity_physical_group_data!(s)
    _organize_pressure_physical_group_data!(s)

    _write_pq_matrices(s)

    cpp_function::Ref{Function} =
    if _is_quadratic(s)
        _cpp_simulate_fem_helmholtz_tetrahedron_quadratic
    else
        _cpp_simulate_fem_helmholtz_tetrahedron_linear
    end
    cpp_function[](
        # freq vector
        s._fi_to_freq,
        s._fi_count,
        # mesh nodes
        s._ni_to_xyz,
        s._ni_count,
        # volume materials
        UInt64.(vec(s._vei_to_ni) .- 1),
        UInt64.(s._vei_to_ivpg .- 1),
        s._vei_count,
        vec(s._ivpg_to_density_values),
        vec(s._ivpg_to_soundspeed_values),
        s._ivpg_count,
        # surface materials
        UInt64.(vec(s._isei_to_ni) .- 1),
        UInt64.(s._isei_to_ispgi .- 1),
        s._isei_count,
        vec(s._ispgi_to_impedance_values),
        s._ispgi_count,
        # volume velocity
        UInt64.(s._vpi_to_ni .- 1),
        vec(s._vpi_to_volvel_values),
        s._vpi_count,
        # surface velocity
        UInt64.(vec(s._vsei_to_ni) .- 1),
        UInt64.(s._vsei_to_ispgv .- 1),
        s._vsei_count,
        vec(s._ispgv_to_velocity_values),
        s._ispgv_count,
        # pressure
        UInt64.(s._pni_to_ni .- 1),
        s._pni_count,
        UInt64.(s._pvi_to_pni_count),
        vec(s._pvi_to_pressure_values),
        s._pvi_count,
        # export
        s._hdf5_file_path
    )
    s._did_run = true
end

function _check_if_did_run(s::Simulation)
    if s._did_run
        error("This simulation has already been run.")
    end
end
