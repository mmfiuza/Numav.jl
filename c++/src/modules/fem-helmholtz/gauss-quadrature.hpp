// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#include "numav/numav.hpp"
#include "common/maths.hpp"

namespace numav {

template<ElementShape S, ElementOrder O> constexpr uint64_t NGP_FORC = [] {
    if constexpr (S==ElementShape::TETRAHEDRON && O==ElementOrder::LINEAR) {
         return 1UL;
    }
    if constexpr (S==ElementShape::TETRAHEDRON && O==ElementOrder::QUADRATIC) {
        return 3UL;
    }
}();

template<ElementShape S, ElementOrder O> constexpr uint64_t NGP_DAMP = [] {
    if constexpr (S==ElementShape::TETRAHEDRON && O==ElementOrder::LINEAR) {
        return 3UL;
    }
    if constexpr (S==ElementShape::TETRAHEDRON && O==ElementOrder::QUADRATIC) {
        return 7UL;
    }
}();

template<ElementShape S, ElementOrder O> constexpr uint64_t NGP_STIF = [] {
    if constexpr (S==ElementShape::TETRAHEDRON && O==ElementOrder::LINEAR) {
        return 1UL;
    }
    if constexpr (S==ElementShape::TETRAHEDRON && O==ElementOrder::QUADRATIC) {
        return 4UL;
    }
}();

template<ElementShape S, ElementOrder O> constexpr uint64_t NGP_MASS = [] {
    if constexpr (S==ElementShape::TETRAHEDRON && O==ElementOrder::LINEAR) {
        return 4UL;
    }
    if constexpr (S==ElementShape::TETRAHEDRON && O==ElementOrder::QUADRATIC) {
        return 15UL;
    }
}();

template<ElementShape S, uint64_t N>
constexpr std::array<std::array<Float,2UL>,N> GAUSS_POINTS_SFC = [] {
    if constexpr (S == ElementShape::TETRAHEDRON) {
        if constexpr (N == 1UL) {
            constexpr Float a = 1_F / 3_F;
            return std::array<std::array<Float,2UL>,N>{{ {a,a} }};
        }
        if constexpr (N == 3UL) {
            constexpr Float a = 1_F / 6_F;
            constexpr Float b = 2_F / 3_F;
            return std::array<std::array<Float,2UL>,N>{{
                {a,a}, {b,a}, {a,b}
            }};
        }
        if constexpr (N == 4UL) {
            constexpr Float a = 1_F / 3_F;
            constexpr Float b = 3_F / 5_F;
            constexpr Float c = 1_F / 5_F;
            return std::array<std::array<Float,2UL>,N>{{
                {a,a}, {b,c}, {c,b}, {c,c}
            }};
        }
        if constexpr (N == 6UL) {
            constexpr Float a  = (6_F - SQRT(15_F)) / 21_F;
            constexpr Float b  = (6_F + SQRT(15_F)) / 21_F;
            constexpr Float c = 1_F - 2_F * a;
            constexpr Float d = 1_F - 2_F * b;
            return std::array<std::array<Float,2UL>,N>{{
                {a,a}, {a,c}, {c,a}, {b,b}, {b,d}, {d,b}
            }};
        }
        if constexpr (N == 7UL) {
            constexpr Float a = 1_F / 3_F;
            constexpr Float b = (6_F + SQRT(15_F)) / 21_F;
            constexpr Float c = (6_F - SQRT(15_F)) / 21_F;
            constexpr Float d = 1_F - 2_F * b;
            constexpr Float e = 1_F - 2_F * c;
            return std::array<std::array<Float,2UL>,N>{{
                {a,a}, {b, b}, {b,d}, {d,b}, {c,c}, {c,e}, {e,c}
            }};
        }
        if constexpr (N == 9UL) {
            constexpr Float a = 0.124949503233232_F; //todo: discover the fracs
            constexpr Float b = 0.437525248383384_F;
            constexpr Float c = 0.797112651860071_F;
            constexpr Float d = 0.165409927389841_F;
            constexpr Float e = 0.037477420750088_F;
            return std::array<std::array<Float,2UL>,N>{{
                {a,b}, {b,a}, {b,b}, {c,d}, {c,e}, {d,c}, {d,e}, {e,c}, {e,d}
            }};
        }
        if constexpr (N == 12UL) {
            constexpr Float a = 0.063089014491502_F; //todo: discover the fracs
            constexpr Float b = 1_F - 2*a;
            constexpr Float c = 0.249286745170910_F;
            constexpr Float d = 1_F - 2*c;
            constexpr Float e = 0.636502499121399_F;
            constexpr Float f = 0.310352451033785_F;
            constexpr Float g = 1_F - e - f;
            return std::array<std::array<Float,2UL>,N>{{
                {b,a}, {a,b}, {a,a}, {d,c}, {c,d}, {c,c},
                {e,f}, {e,g}, {f,e}, {f,g}, {g,e}, {g,f}
            }};
        }
    }
}();

template<ElementShape S, uint64_t N>
constexpr std::array<std::array<Float,3UL>,N> GAUSS_POINTS_VOL = [] {
    if constexpr (S == ElementShape::TETRAHEDRON) {
        if constexpr (N == 1UL) {
            constexpr Float a = 1_F / 4_F;
            return std::array<std::array<Float,3UL>,N>{{ {a, a, a} }};
        }
        if constexpr (N == 4UL) {
            constexpr Float a = (5_F -     SQRT(5_F)) / 20_F;
            constexpr Float b = (5_F + 3_F*SQRT(5_F)) / 20_F;
            return std::array<std::array<Float,3UL>,N>{{
                {a,a,a}, {b,a,a}, {a,b,a}, {a,a,b}
            }};
        }
        if constexpr (N == 5UL) {
            constexpr Float a = 1_F / 4_F;
            constexpr Float b = 1_F / 6_F;
            constexpr Float c = 1_F / 2_F;
            return std::array<std::array<Float,3UL>,N>{{
                {a,a,a},
                {b,b,b}, {c,b,b}, {b,c,b}, {b,b,c}
            }};
        }
        if constexpr (N == 11UL) {
            constexpr Float a =  1_F /  4_F;
            constexpr Float b =  1_F / 14_F;
            constexpr Float c = 11_F / 14_F;
            constexpr Float d = (14_F + SQRT(70_F)) / 56_F;
            constexpr Float e = (14_F - SQRT(70_F)) / 56_F;
            return std::array<std::array<Float,3UL>,N>{{
                {a,a,a},
                {b,b,b}, {c,b,b}, {b,c,b}, {b,b,c},
                {d,d,e}, {d,e,d}, {e,d,d}, {d,e,e}, {e,d,e}, {e,e,d}
            }};
        }
        if constexpr (N == 15UL) {
            constexpr Float a = 1_F /  4_F;
            constexpr Float b = ( 7_F +     SQRT(15_F)) / 34_F;
            constexpr Float c = ( 7_F -     SQRT(15_F)) / 34_F;
            constexpr Float d = (13_F - 3_F*SQRT(15_F)) / 34_F;
            constexpr Float e = (13_F + 3_F*SQRT(15_F)) / 34_F;
            constexpr Float f = ( 5_F -     SQRT(15_F)) / 20_F;
            constexpr Float g = ( 5_F +     SQRT(15_F)) / 20_F;
            return std::array<std::array<Float,3UL>,N>{{
                {a,a,a},
                {b,b,b}, {b,b,d}, {b,d,b}, {d,b,b},
                {c,c,c}, {c,c,e}, {c,e,c}, {e,c,c},
                {f,f,g}, {f,g,f}, {g,f,f}, {f,g,g}, {g,f,g}, {g,g,f}
            }};
        }
    }
}();

template<ElementShape S, uint64_t N>
constexpr std::array<Float,N> GAUSS_WEIGHTS_SFC = [] {
    if constexpr (S == ElementShape::TETRAHEDRON) {
        if constexpr (N == 1UL) {
            constexpr Float a = 1_F / 2_F;
            return std::array<Float,N>({a});
        }
        if constexpr (N == 3UL) {
            constexpr Float a = 1_F / 6_F;
            return std::array<Float,N>({a,a,a});
        }
        if constexpr (N == 4UL) {
            constexpr Float a = -9_F / 32_F;
            constexpr Float b = 25_F / 96_F;
            return std::array<Float,N>({a, b, b, b});
        }
        if constexpr (N == 6UL) {
            constexpr Float a = (8_F*SQRT(15_F) - 15_F) / (96_F*SQRT(15_F));
            constexpr Float b = (8_F*SQRT(15_F) + 15_F) / (96_F*SQRT(15_F));
            return std::array<Float,N>({a, a, a, b, b, b});
        }
        if constexpr (N == 7UL) {
            constexpr Float a = 9_F / 80_F;
            constexpr Float b = (155_F - SQRT(15_F)) / 2400_F;
            constexpr Float c = (155_F + SQRT(15_F)) / 2400_F;
            return std::array<Float,N>({a, c, c, c, b, b, b});
        }
        if constexpr (N == 9UL) {
            constexpr Float a = 0.1029752523804435_F; //todo: discover the fracs
            constexpr Float b = 0.0318457071431115_F;
            return std::array<Float,N>({a, a, a, b, b, b, b, b, b});
        }
        if constexpr (N == 12UL) {
            constexpr Float a = 0.0254224531851035_F; //todo: discover the fracs
            constexpr Float b = 0.0583931378631895_F;
            constexpr Float c = 0.041425537809187_F;
            return std::array<Float,N>({a, a, a, b, b, b, c, c, c, c, c, c});
        }
    }
}();

template<ElementShape S, uint64_t N>
constexpr std::array<Float,N> GAUSS_WEIGHTS_VOL = [] {
    if constexpr (S == ElementShape::TETRAHEDRON) {
        if constexpr (N == 1UL) {
            constexpr Float a = 1_F / 6_F;
            return std::array<Float,N>({a});
        }
        if constexpr (N == 4UL) {
            constexpr Float a = 1_F / 24_F;
            return std::array<Float,N>({a,a,a,a});
        }
        if constexpr (N == 5UL) {
            constexpr Float a = -2_F / 15_F;
            constexpr Float b =  3_F / 40_F;
            return std::array<Float,N>({a,b,b,b,b});
        }
        if constexpr (N == 11UL) {
            constexpr Float a = -74_F /  5625_F;
            constexpr Float b = 343_F / 45000_F;
            constexpr Float c =  28_F /  1125_F;
            return std::array<Float,N>({a,b,b,b,b,c,c,c,c,c,c});
        }
        if constexpr (N == 15UL) {
            constexpr Float a =                        8_F /    405_F;
            constexpr Float b = (2665_F - 14_F*SQRT(15_F)) / 226800_F;
            constexpr Float c = (2665_F + 14_F*SQRT(15_F)) / 226800_F;
            constexpr Float d =                        5_F /    567_F;
            return std::array<Float,N>({a,b,b,b,b,c,c,c,c,d,d,d,d,d,d});
        }
    }
}();

} // namespace numav