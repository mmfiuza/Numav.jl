# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file has the Numav module, that includes all the source code.

module Numav
    
    # automatically create docstrings from files in src/docs
    import ExternalDocstrings
    ExternalDocstrings.@define_docstrings

    include("cpp.jl")
    include("utils.jl")
    include("options.jl")
    include("simulation.jl")
    include("log.jl")
    include("fem-helmholtz/simulation.jl")
    include("fem-helmholtz/materials.jl")
    include("fem-helmholtz/organize-data.jl")
    include("fem-helmholtz/plot.jl")
    include("fem-helmholtz/run.jl")
    include("frequency.jl")
    include("mesh.jl")
    
end
