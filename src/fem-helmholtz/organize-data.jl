# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

function _organize_velocity_region_data!(s::SimulationFemHelmholtz)
    s._vsei_count = 0
    for sei in 1:s._sei_count
        if haskey(s._esr_isrv_bimap, s._sei_to_esr[sei])
            s._vsei_count += 1
        end
    end

    s._vsei_to_ni = Matrix{Int}(undef, _enis_count(s), s._vsei_count)
    s._vsei_to_isrv = Vector{Int}(undef, s._vsei_count)
    vsei = 0
    for sei in 1:s._sei_count
        if haskey(s._esr_isrv_bimap, s._sei_to_esr[sei])
            vsei += 1
            s._vsei_to_ni[:, vsei] = s._sei_to_ni[:, sei]
            s._vsei_to_isrv[vsei] = s._esr_isrv_bimap[s._sei_to_esr[sei]]
        end
    end
    return nothing
end

function _organize_pressure_region_data!(s::SimulationFemHelmholtz)
    s._psei_count = 0
    for sei in 1:s._sei_count
        if haskey(s._esr_isrp_bimap, s._sei_to_esr[sei])
            s._psei_count += 1
        end
    end

    s._psei_to_ni = Matrix{Int}(undef, _enis_count(s), s._psei_count)
    s._psei_to_isrp = Vector{Int}(undef, s._psei_count)
    psei = 0
    for sei in 1:s._sei_count
        if haskey(s._esr_isrp_bimap, s._sei_to_esr[sei])
            psei += 1
            s._psei_to_ni[:, psei] = s._sei_to_ni[:, sei]
            s._psei_to_isrp[psei] = s._esr_isrp_bimap[s._sei_to_esr[sei]]
        end
    end

    # create vector of ni sets
    ni_sets::Vector{Set{Int}} = [Set{Int}() for _ in 1:s._isrp_count]
    for psei in 1:s._psei_count
        ni = s._psei_to_ni[:, psei]
        isrp::Int = s._psei_to_isrp[psei]
        push!(ni_sets[isrp], ni...)
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
            if idx <= s._isrp_count
                push!(funcs, s._isrp_to_pressure_func[idx])
            else
                ppi = idx - s._isrp_count
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

function _organize_volume_region_data!(s::SimulationFemHelmholtz)
    s._vei_to_ivr = Vector{Int}(undef, s._vei_count)
    for vei in 1:s._vei_count
        s._vei_to_ivr[vei] = s._evr_ivr_bimap[s._vei_to_evr[vei]]
    end
    return nothing
end

function _organize_impedance_region_data!(s::SimulationFemHelmholtz)
    s._isei_count = 0
    for sei in 1:s._sei_count
        if haskey(s._esr_isri_bimap, s._sei_to_esr[sei])
            s._isei_count += 1
        end
    end

    s._isei_to_ni = Matrix{Int}(undef, _enis_count(s), s._isei_count)
    s._isei_to_isri = Vector{Int}(undef, s._isei_count)
    isei = 0
    for sei in 1:s._sei_count
        if haskey(s._esr_isri_bimap, s._sei_to_esr[sei])
            isei += 1
            s._isei_to_ni[:, isei] = s._sei_to_ni[:, sei]
            s._isei_to_isri[isei] = s._esr_isri_bimap[s._sei_to_esr[sei]]
        end
    end
    return nothing
end

function _write_pq_matrices(s::SimulationFemHelmholtz)
    # _ivr_to_density_values
    s._ivr_to_density_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._ivr_count, s._fi_count)
    for fi in 1:s._fi_count
        for ivr in 1:s._ivr_count
            s._ivr_to_density_values[ivr, fi] = 
                s._ivr_to_density_func[ivr](s._fi_to_freq[fi])
        end
    end
    
    # _ivr_to_soundspeed_values
    s._ivr_to_soundspeed_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._ivr_count, s._fi_count)
    for ivr in 1:s._ivr_count
        for fi in 1:s._fi_count
            s._ivr_to_soundspeed_values[ivr, fi] =
                s._ivr_to_soundspeed_func[ivr](s._fi_to_freq[fi])
        end
    end

    # _isri_to_impedance_values
    s._isri_to_impedance_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._isri_count, s._fi_count)
    for isri in 1:s._isri_count
        for fi in 1:s._fi_count
            s._isri_to_impedance_values[isri, fi] =
                s._isri_to_impedance_func[isri](s._fi_to_freq[fi])
        end
    end

    # _isrv_to_velocity_values
    s._isrv_to_velocity_values::Matrix{ComplexF64} =
        Matrix{ComplexF64}(undef, s._isrv_count, s._fi_count)
    for isrv in 1:s._isrv_count
        for fi in 1:s._fi_count
            s._isrv_to_velocity_values[isrv, fi] =
                s._isrv_to_velocity_func[isrv](s._fi_to_freq[fi])
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

    return nothing
end
