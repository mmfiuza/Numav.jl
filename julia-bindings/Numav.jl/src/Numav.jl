# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

module Numav

    # automatically create docstrings from files in src/docs
    import ExternalDocstrings
    ExternalDocstrings.@define_docstrings

    include("cpp.jl") # wrapped C++ part
    include("utils.jl") # various things
    include("options.jl") # define every Option singleton

    include("simulation.jl") # general definition of the Simulation type
    include("mesh.jl")
    include("frequency.jl")
    include("fem-helmholtz.jl")
    
end
