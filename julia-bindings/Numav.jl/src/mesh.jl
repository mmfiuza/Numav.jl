# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# export functions

export load_mesh!

function _check_if_simulation_has_mesh(s::Simulation)
    if ( 
        !hasproperty(s, :numerical_method) ||
        !(s.numerical_method isa FemNumericalMethod)
    )
        _throw_simulation_not_applicable()
    end
    return
end

function load_mesh!(
    simulation::Simulation,
    path_to_mesh::AbstractString
)
    _check_if_simulation_has_mesh(simulation)
    _cpp_load_mesh!(simulation._cpp_simulation, String(path_to_mesh))
end
