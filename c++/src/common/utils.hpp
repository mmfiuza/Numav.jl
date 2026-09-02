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

} // namespace numav
