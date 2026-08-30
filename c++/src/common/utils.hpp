// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"

#include <string_view>
#include <charconv>
#include <tuple>
#include <fstream>
#include <memory>

namespace numav {

template<typename T>
std::tuple<T,T> make_ascending_tuple(const T a, const T b) {
    return a<b ? std::make_tuple(a,b) : std::make_tuple(b,a);
}

template<typename T, std::uint64_t... Sizes>
constexpr auto concat_constexpr_arrays(const std::array<T, Sizes>... arrays)
{
    constexpr std::uint64_t total_size = (Sizes + ...);
    std::array<T, total_size> result{};
    std::uint64_t index = 0UL;
    (
        (
            [&]() {
                for (std::uint64_t i = 0UL; i < arrays.size(); ++i) {
                    result[index + i] = arrays[i];
                }
            }(),
            index += arrays.size()
        ), ...
    );
    return result;
}

template<typename T1, typename T2>
std::unique_ptr<T1[]> static_cast_contiguous_data(
    const T2* const data,
    const uint64_t element_count
) {
    std::unique_ptr<T1[]> ptr = std::make_unique<T1[]>(element_count);
    for (uint64_t i = 0UL; i!= element_count; ++i){
        ptr[i] = static_cast<T1>(data[i]);
    }
    return ptr;
}

template< typename T>
std::unique_ptr<uint64_t[]> to_one_based_index(
    const T* const data,
    const uint64_t element_count
) {
    std::unique_ptr<uint64_t[]> ptr = std::make_unique<uint64_t[]>(element_count);
    for (uint64_t i = 0UL; i!= element_count; ++i){
        ptr[i] = data[i] + 1UL;
    }
    return ptr;
}

} // namespace numav
