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

#pragma once

#include "functions/GaussExp.h"
#include "functions/GaussFunc.h"
#include "utils/Printer.h"


namespace mrcpp {

/** @class HeatKernel.
 *
 * @brief Heat kernel in \f$ \mathbb R^D \f$.
 *
 * @details In \f$ \mathbb R^D \f$ the heat kernel has the form
 * \f[
 *  K_t(x)
 *  =
 *   \frac 1{ (4 \pi t)^{D/2} }
 *   \exp
 *   \left(
 *       - \frac{ |x|^2 }{4t}
 *   \right)
 *   , \quad
 *   x \in \mathbb R^D
 *   \text{ and }
 *   t > 0
 *   .
 * \f]
 *
 */
class HeatKernel : public GaussExp<1> {
public:
    HeatKernel()
            : GaussExp<1>() {}

    HeatKernel(double t, int D)
            : GaussExp<1>() {
        constant_coeff = -1.0;

        double expo = 0.25 / t;
        double coef = std::pow(expo / mrcpp::pi, D / 2.0);
        GaussFunc<1> gFunc(expo, coef);
        this->append(gFunc);
    }

    HeatKernel(double t, std::vector<double> &coeffs, int D)
            : GaussExp<1>() {
        initialize(t, coeffs, D);
    }

    void initialize(double t, std::vector<double> &coeffs, int D) {
        constant_coeff = coeffs[0];

        MSG_INFO("Constant: " << constant_coeff);

        for (size_t i = 1; i < coeffs.size(); i++) {
            double exponent = 0.25 / (t * i);
            double coeff = coeffs[i] * std::pow(exponent / mrcpp::pi, D / 2.0);

            MSG_INFO("coeff " << i << ": " << coeffs[i]);

            mrcpp::GaussFunc<1> g(exponent, coeff);

            append(g);
        }
    }

    HeatKernel(double t, int order, int D, std::string &type)
            : GaussExp<1>() {
        if (type == "log") {
            initLogHeatKernel(t, order, 3);
        } else if (type == "ed") {
            initEDHeatKernel(t, order, 3);
        }
    }

    void initEDHeatKernel(double t, int order, int D);
    void initLogHeatKernel(double t, int order, int D);

    double constant_coeff{0.0};
};

} // namespace mrcpp
