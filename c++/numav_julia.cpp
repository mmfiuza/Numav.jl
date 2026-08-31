// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#include "jlcxx/jlcxx.hpp"
#include "numav/numav.hpp"

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.method("_cpp_simulate_fem_helmholtz_tetrahedron_linear",
        &numav::simulate_fem_helmholtz<
            numav::ElementShape::TETRAHEDRON,
            numav::ElementOrder::LINEAR
        >
    );
    mod.method("_cpp_simulate_fem_helmholtz_tetrahedron_quadratic",
        &numav::simulate_fem_helmholtz<
            numav::ElementShape::TETRAHEDRON,
            numav::ElementOrder::QUADRATIC
        >
    );
}
