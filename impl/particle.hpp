#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "./overload.hpp"
#include <array>
#include <cmath>

namespace smd {
class Particle {
public:
    Particle(const double init_x, const double init_y, const double init_z);
    std::array<double,3> coord;
    std::array<double,3> calc_spring(const Particle& other_p, const double k, const double r0);
    std::array<double,3> calc_LennardJones(const Particle& other_p, const double epsilon, const double sigma, const double a);
    std::array<double,3> calc_repulsive(const Particle& other_p, const double a);
    std::array<double,3> calc_soft_repulsive(const Particle& other_p, const double repulsive_d);
private:
    std::array<double,3> calc_force(const Particle& other_p, const double dUdr);
    double spring_dUdr(const double r, const double k, const double r0);
    double LennardJones_dUdr(const double r, const double epsilon, const double sigma, const double a);
    double repulsive_dUdr(const double r, const double a);
    double soft_repulsive_dUdr(const double r, const double repulsive_d);
};

Particle::Particle(const double init_x, const double init_y, const double init_z) {
    coord = {init_x, init_y, init_z};
}

std::array<double,3> Particle::calc_spring(const Particle& other_p, const double k, const double r0) {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = spring_dUdr(r, k, r0);
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_LennardJones(const Particle& other_p, const double epsilon, const double sigma, const double a) {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = LennardJones_dUdr(r, epsilon, sigma, a);
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_repulsive(const Particle& other_p, const double a) {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = repulsive_dUdr(r, a);
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_soft_repulsive(const Particle& other_p, const double repulsive_d) {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = soft_repulsive_dUdr(r, repulsive_d);
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_force(const Particle& other_p, const double dUdr) {
    const auto v_12 = other_p.coord - coord;
    const auto distance = norm(v_12);
    return (-dUdr*(1.0/distance))*v_12;
}

double Particle::spring_dUdr(const double r, const double k, const double r0) {
    return k*(r-r0);
}

double Particle::LennardJones_dUdr(const double r, const double epsilon, const double sigma, const double a) {
    return 24.0*a*epsilon*(2.0*(std::pow(sigma, 12.0)/std::pow(r, 13.0)) - std::pow(sigma, 6.0)/std::pow(r,7.0));
}

double Particle::repulsive_dUdr(const double r, const double a) {
    return -a*(1.0/(r*r));
}

double Particle::soft_repulsive_dUdr(const double r, const double repulsive_d) {
    if (r > repulsive_d) return 0.0;
    return 0.1*(repulsive_d - r);
}
} // smd

#endif