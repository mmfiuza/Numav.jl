# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

function _check_if_it_can_run(s::SimulationFemHelmholtz)
    _check_if_mesh_is_defined(s)
    _check_if_did_run(s)
    if isempty(s._hdf5_file_path)
        error("Result export path not defined." *
              " Call set_result_export_path! to do so.")
    end
    if !s._is_any_source_defined
        error("No sound source was defined." *
              " Call add_sound_source! to do so.")
    end
    if !s._is_freq_defined
        error("Simulation frequency was not defined." *
              " Call set_maximum_frequency! to do so.")
    end
    for evpg in s._existing_evpg
        if !haskey(s._evpg_ivpg_bimap, evpg)
            error("Volume physical group $evpg was not assigned." *
                  " Call add_volume_material! to do so.")
        end
    end
    return nothing
end

function _call_after_every_iteration()
    if !_disable_progress_bar
        _update_progress_bar()
    end
    return nothing
end

function run!(s::SimulationFemHelmholtz)
    _check_if_it_can_run(s)

    _organize_volume_physical_group_data!(s)
    _organize_impedance_physical_group_data!(s)
    _organize_velocity_physical_group_data!(s)
    _organize_pressure_physical_group_data!(s)
    _write_pq_matrices(s)

    _print_simulation_started()
    _print_start_time()
    if !_disable_progress_bar
        bar_tick_count = s._fi_to_freq[1] == 0.0 ? s._fi_count-1 : s._fi_count
        _create_progress_bar(bar_tick_count)
    end

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
        s._hdf5_file_path,
        # other
        CxxWrap.@safe_cfunction(_call_after_every_iteration, Cvoid, ())
    )
    if !_disable_progress_bar
        _finish_progress_bar()
    end
    _print_finish_time()
    s._did_run = true
    return nothing
end
