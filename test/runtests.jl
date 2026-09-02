# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using Test
import Logging
using Numav

@testset "Numav.jl" begin
    Numav._disable_progress_bar = true
    Logging.with_logger(Logging.NullLogger()) do
        include("general.jl")
        include("accuracy.jl")
    end
end
