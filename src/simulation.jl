# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file has code applicable for any Simulation subtype

export
    Simulation,
    create_simulation,
    set_result_export_path!,
    run!

# define general Simulation type
abstract type Simulation end

function create_simulation(;
    numerical_method::Option,
    equation::Option,
    element_shape::Option,
    element_order::Option
)
    if numerical_method === Fem && equation === Helmholtz

        if element_shape === Tetrahedron
            element_shape_type = ElementShape_TETRAHEDRON
        else
            throw(ArgumentError("Invalid `element_shape` option"))
        end

        if element_order === Linear
            element_order_type = ElementOrder_LINEAR
        elseif element_order === Quadratic
            element_order_type = ElementOrder_QUADRATIC
        else
            throw(ArgumentError("Invalid `element_order` option"))
        end

        return SimulationFemHelmholtz{
            element_shape_type,
            element_order_type,
            _enis_count(element_order_type),
            _eniv_count(element_order_type),
        }()
    end
    throw(ArgumentError("Invalid options"))
end

function set_result_export_path!(
    s::Simulation,
    path_to_hdf5_file::AbstractString
)
    _check_if_did_run(s)
    if !isempty(s._hdf5_file_path)
        error("Result export path is already defined.")
    end
    s._hdf5_file_path = path_to_hdf5_file
    return nothing
end

function _check_if_did_run(s::Simulation)
    if s._did_run
        error("This simulation has already been run.")
    end
    return nothing
end
