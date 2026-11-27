# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    add_volume_material!,
    add_surface_material!,
    add_sound_source!

function _check_if_esr_is_valid(s::SimulationFemHelmholtz, esr::Int)
    if !(esr in s._existing_esr)
        error("Region $esr not found in mesh file.")
    end
    if (haskey(s._esr_isri_bimap, esr) ||
        haskey(s._esr_isrv_bimap, esr) ||
        haskey(s._esr_isrp_bimap, esr)
    )
        error("Region $esr already assigned.")
    end
    return nothing
end

function add_volume_material!( 
    s::SimulationFemHelmholtz;
    region::Union{Integer, AbstractVector{<:Integer}},
    density::Fdpq,
    speed_of_sound::Fdpq
)
    if region isa AbstractVector
        for r in region
            add_volume_material!(
                s,
                region = r,
                density = density,
                speed_of_sound = speed_of_sound
            )
        end
        return
    end
    _check_if_mesh_is_defined(s)
    _check_if_did_run(s)
    evr::Int = region
    if !(evr in s._existing_evr)
        error("Region $evr not found in mesh file.")
    end
    if haskey(s._evr_ivr_bimap, evr)
        error("Region $evr already assigned.")
    end
    
    density_func::Function = _fdpq_to_function(density)
    soundspeed_func::Function = _fdpq_to_function(speed_of_sound)
    s._ivr_count += 1
    ivr::Int = s._ivr_count
    s._evr_ivr_bimap[evr] = ivr
    push!(s._ivr_to_density_func, density_func)
    push!(s._ivr_to_soundspeed_func, soundspeed_func)
    return nothing
end

function add_surface_material!( 
    s::SimulationFemHelmholtz;
    region::Union{Integer, AbstractVector{<:Integer}},
    specific_acoustic_impedance::Fdpq
)
    if region isa AbstractVector
        for r in region
            add_surface_material!(
                s,
                region = r,
                specific_acoustic_impedance = specific_acoustic_impedance
            )
        end
        return
    end
    _check_if_mesh_is_defined(s)
    _check_if_did_run(s)
    esr::Int = region
    _check_if_esr_is_valid(s, esr)

    s._isri_count += 1
    isri::Int = s._isri_count
    s._esr_isri_bimap[esr] = isri

    impedance_func::Function = _fdpq_to_function(specific_acoustic_impedance)
    push!(s._isri_to_impedance_func, impedance_func);
    return nothing
end

function add_sound_source!( 
    s::SimulationFemHelmholtz;
    coordinates::Union{
        AbstractVector{<:Real},
        AbstractVector{<:AbstractVector{<:Real}},
        Nothing
    } = nothing,
    region::Union{
        Integer,
        AbstractVector{<:Integer},
        Nothing
    } = nothing,
    volume_velocity::Union{Fdpq, Nothing} = nothing,
    particle_velocity::Union{Fdpq, Nothing} = nothing,
    pressure::Union{Fdpq, Nothing} = nothing
)
    if coordinates isa AbstractVector{<:AbstractVector}
        for c in coordinates
            add_sound_source!(
                s,
                coordinates = c,
                region = region,
                volume_velocity = volume_velocity,
                particle_velocity = particle_velocity,
                pressure = pressure
            )
        end
        return
    end
    if region isa AbstractVector
        for r in region
            add_sound_source!(
                s,
                coordinates = coordinates,
                region = r,
                volume_velocity = volume_velocity,
                particle_velocity = particle_velocity,
                pressure = pressure
            )
        end
        return
    end

    if isnothing(coordinates) && isnothing(region)
        throw(ArgumentError(
            "`coordinates` and `region` not defined"
        ))
    end
    if !isnothing(coordinates) && !isnothing(region)
        throw(ArgumentError(
            "`coordinates` and `region` defined simultaneously"
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
            "x,y,z coordinates do not have 3 components"
        ))
    end
    if !isnothing(volume_velocity) && !isnothing(region)
        throw(ArgumentError(
            "`volume_velocity` and `region` defined simultaneously"
        ))
    end
    if !isnothing(particle_velocity) && !isnothing(coordinates)
        throw(ArgumentError(
            "`particle_velocity` and `coordinates` defined simultaneously"
        ))
    end
    _check_if_mesh_is_defined(s);
    _check_if_did_run(s);
    # TODO: check if the point is outside the mesh

    if !isnothing(coordinates)
        closest_ni::Int = _get_closest_node(s, coordinates)
    elseif !isnothing(region)
        esr::Int = region
        _check_if_esr_is_valid(s, esr);
    end

    # Check if volume_velocity, particle_velocity or pressure was given
    if !isnothing(volume_velocity)
        volvel_func::Function = _fdpq_to_function(volume_velocity)
        push!(s._vpi_to_volvel_func, volvel_func)
        push!(s._vpi_to_ni, closest_ni)
        s._vpi_count += 1
    elseif !isnothing(particle_velocity)
        parvel_func::Function = _fdpq_to_function(particle_velocity)
        push!(s._isrv_to_velocity_func, parvel_func)
        s._isrv_count += 1
        isrv::Int = s._isrv_count
        s._esr_isrv_bimap[esr] = isrv
    elseif !isnothing(pressure)
        pressure_func::Function = _fdpq_to_function(pressure)
        if !isnothing(coordinates)
            push!(s._ppi_to_pressure_func, pressure_func)
            push!(s._ppi_to_ni, closest_ni)
            s._ppi_count += 1
        elseif !isnothing(region)
            push!(s._isrp_to_pressure_func, pressure_func)
            s._isrp_count += 1
            isrp::Int = s._isrp_count
            s._esr_isrp_bimap[esr] = isrp
        end
    end
    s._is_any_source_defined = true;
    return nothing
end
