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

/* Local utility data structures to generate arbitrary order kernels from
 * equation 3.8 in https://doi.org/10.1016/j.jcp.2026.114878
 */

class Polynomial {
public:
    Polynomial()
            : coeffs({0.0}) {}

    Polynomial(std::vector<double> &coeffs)
            : coeffs(coeffs) {}

    size_t get_degree() { return coeffs.size() - 1; }
    size_t get_coeff(size_t i) { return i < coeffs.size() ? coeffs[i] : 0.0; }

    friend Polynomial &operator+=(Polynomial &a, Polynomial &b) {
        while (a.get_degree() < b.get_degree()) a.coeffs.push_back(0.0);

        for (size_t i = 0; i <= b.get_degree(); i++) { a.coeffs[i] += b.coeffs[i]; }

        return a;
    }

    friend Polynomial operator+(Polynomial &a, Polynomial &b) {
        size_t d = std::max(a.get_degree(), b.get_degree());

        std::vector<double> coeffs(d + 1, 0.0);

        for (size_t i = 0; i <= d; i++) { coeffs[i] = a.get_coeff(i) + b.get_coeff(i); }

        return Polynomial(coeffs);
    }

    friend Polynomial &operator*=(Polynomial &a, double b) {
        for (size_t i = 0; i <= a.get_degree(); i++) { a.coeffs[i] *= b; }

        return a;
    }

    friend Polynomial operator*(Polynomial &a, double b) {
        Polynomial new_a = a;

        new_a *= b;

        return new_a;
    }

    friend Polynomial operator*(double b, Polynomial &a) { return a * b; }

    friend Polynomial operator*(Polynomial &a, Polynomial &b) {
        size_t d = a.get_degree() + b.get_degree();

        std::vector<double> coeffs(d + 1, 0.0);

        for (size_t m = 0; m <= a.get_degree(); m++) {
            for (size_t n = 0; n <= b.get_degree(); n++) { coeffs[m + n] += a.coeffs[m] * b.coeffs[n]; }
        }

        return Polynomial(coeffs);
    }

    friend std::ostream &operator<<(std::ostream &os, Polynomial &p) {
        for (size_t i = 0; i <= p.get_degree(); i++) {
            if (i == 0) {
                os << p.coeffs[0];
            } else {
                if (p.coeffs[i] < 0.0) {
                    os << " - " << -p.coeffs[i];
                } else {
                    os << " + " << p.coeffs[i];
                }

                os << " x^" << i;
            }
        }

        return os;
    }

    std::vector<double> coeffs;
};

class HeatKernelData {
public:
    Polynomial &get_coeffs(size_t order) {
        if (data.size() >= order && data[order - 1].size() > 0) { return data[order - 1][0]; }

        if (order == 1) {
            if (data.size() == 0) { data.push_back(std::vector<Polynomial>()); }

            std::vector<double> coeffs{-1.0, 1.0};

            data[0].push_back(Polynomial(coeffs));

            return data[0][0];
        } else if (order == 2) {
            while (data.size() < 2) { data.push_back(std::vector<Polynomial>()); }

            std::vector<double> coeffs{-1.5, 2.0, -0.5};

            data[1].push_back(Polynomial(coeffs));

            return data[1][0];
        }

        Polynomial new_coeffs = get_coeffs(2);

        double two_pow_m = 2.0;
        double m_fac = 2.0;

        for (size_t m = 3; m <= order; m++) {
            two_pow_m *= 2.0;
            m_fac *= m;

            double prefac = (two_pow_m - 2.0) / m_fac;

            Polynomial tmp = get_coeffs(order - m + 1, m);
            tmp *= prefac;

            new_coeffs += tmp;
        }

        while (data.size() < order) { data.push_back(std::vector<Polynomial>()); }

        data[order - 1].push_back(new_coeffs);

        return data[order - 1][0];
    }

    Polynomial &get_coeffs(size_t order, size_t power) {
        if (power == 1) { return get_coeffs(order); }

        if (data.size() >= order && data[order - 1].size() >= power) { return data[order - 1][power - 1]; }

        Polynomial new_poly = get_coeffs(order, power - 1) * get_coeffs(order);

        data[order - 1].push_back(new_poly);

        return data[order - 1][power - 1];
    }

    std::vector<std::vector<Polynomial>> data;
};

HeatKernel EDHeatKernel(double t, int order, int D) {
    // Assuming generating the coefficient data is cheap and that
    // we will not be constructing many such kernels. If this becomes
    // a problem, we can think about caching this data variable.
    HeatKernelData data;

    Polynomial &coeffs = data.get_coeffs(order);

    return mrcpp::HeatKernel(t, coeffs.coeffs, D);
}

} // namespace mrcpp
