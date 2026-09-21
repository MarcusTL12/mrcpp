/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "HeatKernel.h"

namespace mrcpp {

static double harmonic_number(size_t n) {
    double acc = 0.0;

    // Summing backwards to reduce floating point error (adding small numbers first)
    // This is not super important, but for n = 1000000000 (1e9), the difference of
    // summing forward and backward is ~1e-12
    for (size_t i = n; i > 0; i--) { acc += 1.0 / i; }

    return acc;
}

HeatKernel LogHeatKernel(double t, int order, int D) {
    std::vector<double> coeffs;

    coeffs.push_back(-harmonic_number(order));

    double binom = 1.0;

    for (size_t m = 1; m <= order; m++) {
        binom *= static_cast<double>(order - m + 1) / static_cast<double>(m);
        double c = binom / m;
        coeffs.push_back(m % 2 == 0 ? -c : c);
    }

    return mrcpp::HeatKernel(t, coeffs, D);
}

} // namespace mrcpp
