# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file has the Numav module, that includes all the source code.

module Numav

    # automatically create docstrings from files in src/docs
    import ExternalDocstrings
    ExternalDocstrings.@define_docstrings

    include("cpp.jl") # wrapped C++ part
    include("utils.jl") # various things
    include("options.jl") # define every Option singleton
    include("simulation.jl") # definition of the Simulation type and subtypes

    # include functions
    include("functions/any-simulation.jl")
    include("functions/mesh.jl")
    include("functions/frequency.jl")
    include("functions/fem-helmholtz.jl")
    
end
