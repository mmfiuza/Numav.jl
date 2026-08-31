# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    add_volume_material!,
    add_surface_material!,
    add_sound_source!

function _check_if_espg_is_valid(s::SimulationFemHelmholtz, espg::Int)
    if !(espg in s._existing_espg)
        error("Physical group $espg not found in mesh file.")
    end
    if (haskey(s._espg_ispgi_bimap, espg) ||
        haskey(s._espg_ispgv_bimap, espg) ||
        haskey(s._espg_ispgp_bimap, espg)
    )
        error("Physical group $espg already assigned.")
    end
    return nothing
end

function add_volume_material!( 
    s::SimulationFemHelmholtz;
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    density::Fdpq,
    speed_of_sound::Fdpq
)
    if physical_group isa AbstractVector
        for pg in physical_group
            add_volume_material!(
                s,
                physical_group = pg,
                density = density,
                speed_of_sound = speed_of_sound
            )
        end
        return
    end
    _check_if_mesh_is_defined(s)
    _check_if_did_run(s)
    evpg::Int = physical_group
    if !(evpg in s._existing_evpg)
        error("Physical group $evpg not found in mesh file.")
    end
    if haskey(s._evpg_ivpg_bimap, evpg)
        error("Physical group $evpg already assigned.")
    end
    
    density_func::Function = _fdpq_to_function(density)
    soundspeed_func::Function = _fdpq_to_function(speed_of_sound)
    s._ivpg_count += 1
    ivpg::Int = s._ivpg_count
    s._evpg_ivpg_bimap[evpg] = ivpg
    push!(s._ivpg_to_density_func, density_func)
    push!(s._ivpg_to_soundspeed_func, soundspeed_func)
    return nothing
end

function add_surface_material!( 
    s::SimulationFemHelmholtz;
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    specific_acoustic_impedance::Fdpq
)
    if physical_group isa AbstractVector
        for pg in physical_group
            add_surface_material!(
                s,
                physical_group = pg,
                specific_acoustic_impedance = specific_acoustic_impedance
            )
        end
        return
    end
    _check_if_mesh_is_defined(s)
    _check_if_did_run(s)
    espg::Int = physical_group
    _check_if_espg_is_valid(s, espg)

    s._ispgi_count += 1
    ispgi::Int = s._ispgi_count
    s._espg_ispgi_bimap[espg] = ispgi

    impedance_func::Function = _fdpq_to_function(specific_acoustic_impedance)
    push!(s._ispgi_to_impedance_func, impedance_func);
    return nothing
end

function add_sound_source!( 
    s::SimulationFemHelmholtz;
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
    if coordinates isa AbstractVector{<:AbstractVector}
        for c in coordinates
            add_sound_source!(
                s,
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
                s,
                coordinates = coordinates,
                physical_group = pg,
                volume_velocity = volume_velocity,
                particle_velocity = particle_velocity,
                pressure = pressure
            )
        end
        return
    end

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
            "x,y,z coordinates do not have 3 components"
        ))
    end
    if !isnothing(volume_velocity) && !isnothing(physical_group)
        throw(ArgumentError(
            "`volume_velocity` and `physical_group` defined simultaneously"
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
    elseif !isnothing(physical_group)
        espg::Int = physical_group
        _check_if_espg_is_valid(s, espg);
    end

    # Check if volume_velocity, particle_velocity or pressure was given
    if !isnothing(volume_velocity)
        volvel_func::Function = _fdpq_to_function(volume_velocity)
        push!(s._vpi_to_volvel_func, volvel_func)
        push!(s._vpi_to_ni, closest_ni)
        s._vpi_count += 1
    elseif !isnothing(particle_velocity)
        parvel_func::Function = _fdpq_to_function(particle_velocity)
        push!(s._ispgv_to_velocity_func, parvel_func)
        s._ispgv_count += 1
        ispgv::Int = s._ispgv_count
        s._espg_ispgv_bimap[espg] = ispgv
    elseif !isnothing(pressure)
        pressure_func::Function = _fdpq_to_function(pressure)
        if !isnothing(coordinates)
            push!(s._ppi_to_pressure_func, pressure_func)
            push!(s._ppi_to_ni, closest_ni)
            s._ppi_count += 1
        elseif !isnothing(physical_group)
            push!(s._ispgp_to_pressure_func, pressure_func)
            s._ispgp_count += 1
            ispgp::Int = s._ispgp_count
            s._espg_ispgp_bimap[espg] = ispgp
        end
    end
    s._is_any_source_defined = true;
    return nothing
end
