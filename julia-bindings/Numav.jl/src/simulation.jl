# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    Simulation,
    create_simulation,
    set_result_export_path!,
    run!

abstract type Simulation end

abstract type NumericalMethod end
struct FemNumericalMethod <: NumericalMethod end

abstract type Equation end
struct HelmholtzEquation <: Equation end

abstract type ElementShape end
struct TetrahedronElementShape <: ElementShape end

abstract type ElementOrder end
struct LinearElementOrder <: ElementOrder end
struct QuadraticElementOrder <: ElementOrder end

function _enis_count(element_order::ElementOrder)
    if element_order isa LinearElementOrder
        return 3
    elseif element_order isa QuadraticElementOrder
        return 4
    end
end

function _eniv_count(element_order::ElementOrder)
    if element_order isa LinearElementOrder
        return 4
    elseif element_order isa QuadraticElementOrder
        return 10
    end
end

function create_simulation(;
    numerical_method::Option,
    equation::Option,
    element_shape::Option,
    element_order::Option
)
    if numerical_method === Fem && equation === Helmholtz
        cpp_args = (_cpp_NumericalMethod_fem, _cpp_Equation_helmholtz)

        if element_shape === Tetrahedron
            element_shape_type = TetrahedronElementShape
            cpp_args = (cpp_args..., _cpp_ElementShape_tetrahedron)
        else
            throw(ArgumentError("Invalid `element_shape` option"))
        end

        if element_order === Linear
            element_order_type = LinearElementOrder
            cpp_args = (cpp_args..., _cpp_ElementOrder_linear)
        elseif element_order === Quadratic
            element_order_type = QuadraticElementOrder
            cpp_args = (cpp_args..., _cpp_ElementOrder_quadratic)
        else
            throw(ArgumentError("Invalid `element_order` option"))
        end

        return SimulationFemHelmholtz{
                element_shape_type,
                element_order_type,
                _enis_count(element_order_type()),
                _eniv_count(element_order_type()),
            }( _cpp_simulation = _cpp_Simulation{cpp_args...}() )
    end
    throw(ArgumentError("Invalid options"))
end

function set_result_export_path!(
    simulation::Simulation,
    path_to_hdf5_file::AbstractString
)
    _cpp_set_result_export_path!(
        simulation._cpp_simulation,
        String(path_to_hdf5_file)
    )
end

function run!(
    simulation::Simulation
)
    _cpp_run!(simulation._cpp_simulation)
end
