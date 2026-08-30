# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This files has functions applicable for Simulation subtypes the run iterating
# on frequency steps

export
    set_frequency!,
    get_frequency_vector

const SimulationWithFreqDomain = Union{SimulationFemHelmholtz}

function set_frequency!( 
    s::SimulationWithFreqDomain;
    max::Union{Real, Nothing} = nothing,
    min::Union{Real, Nothing} = nothing,
    length::Union{Integer, Nothing} = nothing,
    sampling_density::Union{Option, Nothing} = nothing,
    step::Union{Real, Nothing} = nothing,
    vector::Union{AbstractVector{<:Real}, Nothing} = nothing
)
    if isnothing(max) && isnothing(vector)
        throw(ArgumentError("Neither `max` nor `vector` passed"))
    end
    if !isnothing(vector)
        if (
            !isnothing(max) ||
            !isnothing(min) ||
            !isnothing(length) ||
            !isnothing(sampling_density) ||
            !isnothing(step)
        )
            throw(ArgumentError(
                "`vector` and other parameter(s) passed simultaneously"
            ))
        end
    end
    if !isnothing(step) && !isnothing(length)
        throw(ArgumentError(
            "`step` and `length` passed simultaneously"
        ))
    end
    if !isnothing(step) && !isnothing(sampling_density)
        throw(ArgumentError(
            "`step` and `sampling_density` passed simultaneously"
        ))
    end
    if s._is_freq_defined
        error("Frequency is already defined.")
    end
    if !isnothing(max)
        min = something(min, 0)
        if !isnothing(step)
            vector = range(min, max, step=step)
        else
            length = something(length, 4096)
            sampling_density = something(sampling_density, Quadratic)
            if sampling_density === Constant
                vector = range(min, max, length=length)
            elseif sampling_density === Quadratic
                vector = _cubic_range(min, max, length)
            else
                throw(ArgumentError("Invalid `sampling_density` option"))
            end
        end
    end
    s._fi_to_freq = Float64.(Vector(vector))
    s._fi_count = Main.length(s._fi_to_freq)
    s._is_freq_defined = true
    return nothing
end

function get_frequency_vector(
    s::SimulationWithFreqDomain
)
    return copy(s._fi_to_freq)
end
