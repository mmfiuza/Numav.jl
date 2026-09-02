# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using HDF5

const HDF5_SIGNATURE = "numav_result_hdf5_0.3.0"
const HDF5_NUMERICAL_METHOD = "finite_element_method"
const HDF5_EQUATION = "helmholtz"

function _hdf5_element_shape_str(shape::ElementShape)
    shape === ElementShape_TETRAHEDRON && return "tetrahedron"
    error("unhandled ElementShape $shape")
end

function _hdf5_element_order_str(order::ElementOrder)
    order === ElementOrder_LINEAR    && return "linear"
    order === ElementOrder_QUADRATIC && return "quadratic"
    error("unhandled ElementOrder $order")
end

function _begin_hdf5_file(s::SimulationFemHelmholtz)
    file = HDF5.h5open(s._hdf5_file_path, "w")
    HDF5.attributes(file)["Conventions"] = HDF5_SIGNATURE

    sim_type_grp = HDF5.create_group(file, "simulation_type")
    HDF5.attributes(sim_type_grp)["numerical_method"] = HDF5_NUMERICAL_METHOD
    HDF5.attributes(sim_type_grp)["equation"] = HDF5_EQUATION
    HDF5.attributes(sim_type_grp)["element_shape"] =
        _hdf5_element_shape_str(typeof(s).parameters[1])
    HDF5.attributes(sim_type_grp)["element_order"] =
        _hdf5_element_order_str(typeof(s).parameters[2])

    solution_grp = HDF5.create_group(file, "results")

    # Julia shape: (ni_count, fi_count)  -> disk/h5dump shape: (fi_count, ni_count)
    pressure_dataset = HDF5.create_dataset(
        solution_grp, "pressure", ComplexF64, (s._ni_count, s._fi_count);
        chunk = (s._ni_count, 1),
        compress = 6,
    )
    HDF5.attributes(pressure_dataset)["units"] = "Pa"

    HDF5.API.h5ds_set_label(pressure_dataset.id, 0,
        Vector{UInt8}(codeunits("frequency_index\0"))
    )
    HDF5.API.h5ds_set_label(pressure_dataset.id, 1,
        Vector{UInt8}(codeunits("node_index\0"))
    )
    return file, pressure_dataset
end

function _write_simulation_inputs_to_hdf5_file(
    s::SimulationFemHelmholtz,
    file::HDF5.File
)
    inputs_grp = HDF5.create_group(file, "inputs")

    # simulated frequencies
    freq_data = Vector{Float64}(s._fi_to_freq)
    dataset_freq = HDF5.create_dataset(inputs_grp, "simulated_frequencies",
        Float64, (length(freq_data),))
    write(dataset_freq, freq_data)
    HDF5.attributes(dataset_freq)["units"] = "Hz"
    HDF5.API.h5ds_set_label(dataset_freq.id, 0,
        Vector{UInt8}(codeunits("frequency_index\0")))
    HDF5.API.h5ds_set_scale(dataset_freq.id, "simulated_frequencies")

    # mesh
    mesh_grp = HDF5.create_group(inputs_grp, "mesh")

    # nodes
    nodes_data = Matrix{Float64}(s._ni_to_xyz)
    dataset_nodes = HDF5.create_dataset(mesh_grp, "nodes", Float64, size(nodes_data))
    write(dataset_nodes, nodes_data)
    HDF5.attributes(dataset_nodes)["units"] = "m"
    HDF5.API.h5ds_set_label(dataset_nodes.id, 0,
        Vector{UInt8}(codeunits("node_index\0")))
    HDF5.API.h5ds_set_label(dataset_nodes.id, 1,
        Vector{UInt8}(codeunits("coordinates\0")))

    # volume_elements
    elem_data = Matrix{UInt64}(s._vei_to_ni)
    dataset_elem = HDF5.create_dataset(mesh_grp, "volume_elements",
        UInt64, size(elem_data))
    write(dataset_elem, elem_data)
    HDF5.attributes(dataset_elem)["units"] = "1"
    HDF5.API.h5ds_set_label(dataset_elem.id, 0,
        Vector{UInt8}(codeunits("volume_element_index\0")))
    HDF5.API.h5ds_set_label(dataset_elem.id, 1,
        Vector{UInt8}(codeunits("elemental_node_index\0")))
    HDF5.attributes(dataset_elem)["node_index_base"] = UInt64(1)

    return nothing
end
