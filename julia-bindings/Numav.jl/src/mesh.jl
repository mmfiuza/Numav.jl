# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# export functions

export load_mesh!

const SimulationWithMesh = Union{SimulationFemHelmholtz}

function load_mesh!(
    simulation::SimulationWithMesh,
    path_to_mesh::AbstractString
)
    _cpp_load_mesh!(simulation._cpp_simulation, String(path_to_mesh))
end
