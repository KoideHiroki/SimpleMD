#ifndef OVERLOAD_HPP
#define OVERLOAD_HPP

#include <array>
#include <cmath>

namespace smd {
    std::array<double,3> operator*(const double a, const std::array<double,3>& v) {
        return {a*v[0], a*v[1], a*v[2]};
    }

    std::array<double,3> operator-(const std::array<double,3>& v1, const std::array<double,3>& v2) {
        return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
    }

    double norm(const std::array<double,3>& v) {
        return std::sqrt(
            v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
        );
    }

    std::array<double,3> operator+=(const std::array<double,3>& lhs, const std::array<double,3>& rhs) {
        return {lhs[0]+rhs[0], lhs[1]+rhs[1], lhs[2]+rhs[2]};
    }

    std::array<double,3> operator/(const std::array<double,3>& v, const double d) {
        return {v[0]/d, v[1]/d, v[2]/d};
    }

    std::array<double,3> operator+(const std::array<double,3>& lhs, const std::array<double,3>& rhs) {
        return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
    }
}


#endif