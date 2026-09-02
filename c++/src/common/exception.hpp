// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace numav {

    template<typename... T>
    void error(T&&... args) {
        std::ostringstream oss;
        (oss << ... << std::forward<T>(args));
        std::string msg = oss.str();
        std::cerr << msg << "\n";
        throw std::runtime_error("numav error");
    }

} // namespace numav
