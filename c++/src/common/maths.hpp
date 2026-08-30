// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "modules/fem-helmholtz/fem-helmholtz.hpp"
#include "common/utils.hpp"

#include "SafePtr.hpp"

#include <tuple>
#include <cmath>

namespace numav {

Float get_triangle_area(
    const std::array<std::array<Float,3UL>,3UL> coords
);

Float get_tetrahedron_volume(
    const std::array<std::array<Float,3UL>,4UL> coords
);

template<uint64_t N> constexpr uint64_t FACTORIAL = [] {
    uint64_t result = 1UL;
    for (uint64_t i = 2UL; i <= N; ++i) {
        result *= i;
    }
    return result;
}();

template<uint64_t K, typename T>
constexpr T power(T N) {
    if constexpr (K == 0) {
        return T(1);
    }
    else if constexpr (K == 1) {
        return N;
    }
    else {
        return N * power<K-1>(N);
    }
}

constexpr Float SQRT(const Float X) {
    Float low_bound = 0;
    Float high_bound = (X > 1) ? X : 1;
    for (uint64_t i = 0UL; i != 100UL; ++i)
    {
        Float mid = (low_bound + high_bound) / 2_F;
        if (mid * mid > X) {
            high_bound = mid;
        }
        else {
            low_bound = mid;
        }
    }
    return low_bound;
};

template<auto X, uint64_t N>
constexpr auto POWER = [] {
    using T = decltype(X);
    T result = T(1);
    for (uint64_t i = 0UL; i != N; ++i) {
        result *= X;
    }
    return result;
}();

template<uint64_t N, uint64_t K>
constexpr uint64_t PERMUTATION_REP_SIZE = [] {
    // combination with repetition
    return POWER<N,K>;
}();

template<uint64_t N>
constexpr std::array<std::array<uint64_t,2UL>, PERMUTATION_REP_SIZE<N,2UL>> 
PERMUTATION_REP = [] { // todo: generalize for K!=2 (maybe)
    constexpr uint64_t K = 2UL;
    // combination with repetition in upper column major order
    std::array<std::array<uint64_t,K>, PERMUTATION_REP_SIZE<N,K>> result;
    auto it_result = result.begin();
    for (uint64_t j = 0UL; j != N; ++j) {
        for (uint64_t i = 0UL; i != N; ++i) {
            *it_result = {i, j};
            ++it_result;
        }
    }
    return result;
}();

template<uint64_t N, uint64_t K> constexpr uint64_t COMBINATION_REP_SIZE = [] {
    // combination with repetition
    return FACTORIAL<N+K-1UL> / (FACTORIAL<K> * FACTORIAL<N-1UL>);
}();

template<uint64_t N>
constexpr std::array<std::array<uint64_t,2UL>, COMBINATION_REP_SIZE<N,2UL>> 
COMBINATION_REP = [] { // todo: generalize for K!=2 (maybe)
    constexpr uint64_t K = 2UL;
    // combination with repetition in upper column major order
    std::array<std::array<uint64_t,K>, COMBINATION_REP_SIZE<N,K>> result;
    auto it_result = result.begin();
    for (uint64_t j = 0UL; j != N; ++j) {
        for (uint64_t i = 0UL; i != j+1UL; ++i) {
            *it_result = {i, j};
            ++it_result;
        }
    }
    return result;
}();

} // namespace numav
