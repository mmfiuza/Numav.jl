# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

@testset "accuracy" begin

import DelimitedFiles
import HDF5

function read_expected_result(file_path::String)
    content = read(file_path, String)
    content = replace(content, "i" => "im")
    fixed_lines = String[]
    for line in split(content, '\n')
        isempty(strip(line)) && continue
        fields = split(line, ',')
        # Inlined normalize_complex
        fixed_fields = map(fields) do f
            s = strip(f)
            return occursin("im", s) ? s : s * "+0.0im"
        end
        push!(fixed_lines, join(fixed_fields, ','))
    end
    fixed_content = join(fixed_lines, '\n')
    data = DelimitedFiles.readdlm(IOBuffer(fixed_content), ',', ComplexF64)
    p = data[:, 4:end]
    xyz = data[:, 1:3]
    return p, xyz
end

function test(linear_or_quadratic::String)
    if linear_or_quadratic == "linear"
        ord = Linear
    elseif linear_or_quadratic == "quadratic"
        ord = Quadratic
    end
    f(x) = x * (1 + 1im)
    # simulate
    s = create_simulation(
        numerical_method = Fem,
        equation = Helmholtz,
        element_shape = Tetrahedron,
        element_order = ord
    )
    set_frequency!(s, min=10, max=100, step=1)
    load_mesh!(s, "test1.bdf")
    add_volume_material!(s, physical_group=1, density=1.2, speed_of_sound=343)
    add_volume_material!(s, physical_group=2, density=1.0, speed_of_sound=100)
    add_volume_material!(s, physical_group=3, density=f, speed_of_sound=f)
    add_surface_material!(s, physical_group=4, specific_acoustic_impedance=400)
    add_surface_material!(s, physical_group=3, specific_acoustic_impedance=f)
    add_sound_source!(s, physical_group=6, pressure=8)
    add_sound_source!(s, physical_group=5, pressure=f)
    add_sound_source!(s, coordinates=[5, 0, 0], pressure=2)
    add_sound_source!(s, coordinates=[5, 4, 0], pressure=f)
    add_sound_source!(s, coordinates=[0, 0, 0], volume_velocity=10)
    add_sound_source!(s, coordinates=[0, 4, 0], volume_velocity=f)
    add_sound_source!(s, physical_group=1, particle_velocity=5)
    add_sound_source!(s, physical_group=2, particle_velocity=f)
    set_result_export_path!(s, "result-actual-"*linear_or_quadratic*".h5")
    run!(s)

    # read simulated result
    h5_actual = HDF5.h5open("result-actual-"*linear_or_quadratic*".h5")
    ni_count = 42
    p_actual = h5_actual["/results/pressure"][1:ni_count, :]
    xyz_actual = h5_actual["/inputs/mesh/nodes"][:, 1:ni_count]

    # reorder nodes
    perm = sortperm(1:size(xyz_actual, 2), by = i -> 
        (xyz_actual[1, i], xyz_actual[2, i], xyz_actual[3, i])
    )
    xyz_actual = xyz_actual[:, perm]
    p_actual = p_actual[perm, :]

    # read expected result
    p_expected, xyz_expected = read_expected_result(
        "p-expected-"*linear_or_quadratic*".txt"
    )

    spl_actual = 20*log10.(abs.(p_actual)/sqrt(2)/20e-6);
    spl_expected = 20*log10.(abs.(p_expected)/sqrt(2)/20e-6);
    spl_diff = maximum(abs.(spl_expected - spl_actual))

    phase_actual = angle.(p_actual);
    phase_expected = angle.(p_expected);
    phase_diff = maximum(abs.(phase_expected - phase_actual))
    
    return spl_diff, phase_diff
end

spl_diff_linear, phase_diff_linear = test("linear")
@test spl_diff_linear < 0.01 && phase_diff_linear < 0.01

spl_diff_quadratic, phase_diff_quadratic = test("quadratic")
@test spl_diff_quadratic < 0.01 && phase_diff_quadratic < 0.01

end # @testset
