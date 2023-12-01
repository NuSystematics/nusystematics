#pragma once

#include <array>

namespace nusyst {
static std::array<double, 6> const Gauss2DParams_CV = {{
    10.5798, 0.254032, 0.50834, 0.0571035, 0.129051, 0.875287}};

static std::array<double, 6> const Gauss2DParams_NNOnly = {{
    17.0344, 0.289916, 0.532062, 0.0746849, 0.137321, 0.836689}};

static std::array<double, 6> const Gauss2DParams_npOnly = {{
    8.58724, 0.23626, 0.502603, 0.072291, 0.154832, 0.789796}};

static std::array<double, 6> const Gauss2DParams_1p1hOnly = {{
    5.38719, 0.213611, 0.396552, 0.0496312, 0.125062, 0.806659}};

static std::array<double, 10> const RPAPolyQ2_CV = {{
    0.578908, 1.36809,     -1.27758,    0.57941,     -0.146737,
    0.021431, -0.00170815, 5.57246e-05, 8.71718e-07, -7.68945e-08}};

static std::array<double, 10> const RPAPolyQ2_Plus1 = {{
    0.578908, 1.36809,     -1.27758,    0.57941,     -0.146737,
    0.021431, -0.00170815, 5.57246e-05, 8.71718e-07, -7.68945e-08}};

static std::array<double, 10> const RPAPolyQ2_Minus1 = {{
    0.578908, 1.36809,     -1.27758,    0.57941,     -0.146737,
    0.021431, -0.00170815, 5.57246e-05, 8.71718e-07, -7.68945e-08}};
} // namespace nusyst
