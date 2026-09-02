// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"

#include <stdexcept>
#include <format>
#include <string>

namespace numav {

    template<typename... T>
    void error(std::format_string<T...> msg_fmt, T&&... args) {
        std::string msg = std::format(msg_fmt, std::forward<T>(args)...);
        std::cerr << msg << "\n";
        throw std::runtime_error("numav error");
    }

} // namespace numav
