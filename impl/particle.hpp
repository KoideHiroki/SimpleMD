#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "./overload.hpp"
#include "./constants.hpp"
#include <array>
#include <cmath>

namespace smd {
class Particle {
public:
    Particle(const double init_x, const double init_y, const double init_z, const double init_weight, const double init_random_x, const double init_random_y, const double init_random_z);
    std::array<double,3> coord;
    std::array<double,3> velo;
    double weight;
    std::array<double,3> random_memory;
    std::array<double,3> calc_spring(const Particle& other_p, const double k, const double r0) const;
    std::array<double,3> calc_LennardJones(const Particle& other_p, const double epsilon, const double sigma) const;
    std::array<double,3> calc_repulsive(const Particle& other_p, const double a) const;
    std::array<double,3> calc_soft_repulsive(const Particle& other_p, const double repulsive_d) const;
    std::array<double,3> calc_excluded(const Particle& other_p, const double excluded_d) const;
    std::array<double,3> calc_sphere() const;
private:
    std::array<double,3> calc_force(const Particle& other_p, const double dUdr) const;
    double spring_dUdr(const double r, const double k, const double r0) const;
    double LennardJones_dUdr(const double r, const double epsilon, const double sigma) const;
    double repulsive_dUdr(const double r, const double a) const;
    double soft_repulsive_dUdr(const double r, const double repulsive_d) const;
    double excluded_dUdr(const double r, const double excluded_d) const;
    double sphere_dUdr() const;

    double debug = 100.0;
};

Particle::Particle(const double init_x, const double init_y, const double init_z, const double init_weight, const double init_random_x, const double init_random_y, const double init_random_z) {
    coord = {init_x, init_y, init_z};
    velo = {0.0, 0.0, 0.0};
    weight = init_weight;
    random_memory = {init_random_x, init_random_y, init_random_z};
}

std::array<double,3> Particle::calc_spring(const Particle& other_p, const double k, const double r0) const {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = spring_dUdr(r, k, r0);
    //const auto force = calc_force(other_p, dUdr);
    //if (norm(force) > debug) {
    //    std::cerr << r << std::endl;
    //    std::cerr << "too big force in spring" << std::endl; 
    //    std::exit(1);
    //}
    //return force;
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_LennardJones(const Particle& other_p, const double epsilon, const double sigma) const {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = LennardJones_dUdr(r, epsilon, sigma);
    //const auto force = calc_force(other_p, dUdr);
    //if (norm(force) > debug) {
    //    std::cerr << r << std::endl;
    //    std::cerr << "too big force in LJ" << std::endl; 
    //    std::exit(1);
    //}
    //return force;
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_repulsive(const Particle& other_p, const double a) const {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = repulsive_dUdr(r, a);
    //const auto force = calc_force(other_p, dUdr);
    //if (norm(force) > debug) {
    //    std::cerr << r << std::endl;
    //    std::cerr << "too big force in rep" << std::endl; 
    //    std::exit(1);
    //}
    //return force;
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_soft_repulsive(const Particle& other_p, const double repulsive_d) const {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = soft_repulsive_dUdr(r, repulsive_d);
    //const auto force = calc_force(other_p, dUdr);
    //if (norm(force) > debug) {
    //    std::cerr << r << std::endl;
    //    std::cerr << "too big force in soft rep" << std::endl; 
    //    std::exit(1);
    //}
    //return force;
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_excluded(const Particle& other_p, const double excluded_d) const {
    const auto r = norm(coord - other_p.coord);
    const auto dUdr = excluded_dUdr(r, excluded_d);
    //const auto force = calc_force(other_p, dUdr);
    //if (norm(force) > debug) {
    //    std::cerr << r << std::endl;
    //    std::cerr << "too big force in excluded" << std::endl;
    //    std::exit(1);
    //}
    //return force;
    return calc_force(other_p, dUdr);
}

std::array<double,3> Particle::calc_sphere() const {
    const auto v = coord;
    const auto n = norm(v);
    const auto dUdr = sphere_dUdr();
    const auto force = -dUdr*(v/n);
    return force;
}

std::array<double,3> Particle::calc_force(const Particle& other_p, const double dUdr) const {
    const auto v_12 = coord - other_p.coord;
    const auto distance = norm(v_12);
    auto force = (-dUdr*(1.0/(distance+1e-6)))*v_12;

    const double f_n = norm(force);
    if (f_n > particle::FMAX) force = (particle::FMAX / f_n) * force;
    return force;
}

double Particle::spring_dUdr(const double r, const double k, const double r0) const {
    return k*(r-r0);
}

double Particle::LennardJones_dUdr(const double r, const double epsilon, const double sigma) const {
    //return 24.0*epsilon*(-2.0*(std::pow(sigma, 12.0)/(std::pow(r, 13.0)+1e-6)) + std::pow(sigma, 6.0)/(std::pow(r, 7.0)+1e-6));
    double inv_r  = 1.0 / (r+1e-6);
    double sr     = sigma * inv_r;
    double sr2    = sr * sr;
    double sr6    = sr2 * sr2 * sr2;
    double sr12   = sr6 * sr6;
    double dUdr   = 24.0 * epsilon * (-2*sr12 + sr6) * inv_r;
    return dUdr;
}

double Particle::repulsive_dUdr(const double r, const double a) const {
    if (r > particle::REPULSIVE_D) return 0.0;
    return -a*(1.0/(r*r+1e-6));
}

double Particle::soft_repulsive_dUdr(const double r, const double repulsive_d) const {
    if (r > repulsive_d) return 0.0;
    return -particle::SOFT_REPULSIVE_A*repulsive_d;
}

double Particle::excluded_dUdr(const double r, const double excluded_d) const {
    if (r > excluded_d) return 0.0;
    return -(1.0+std::pow(std::tan(((3.1415926535/2.0)/excluded_d)*(r-excluded_d)), 2.0));
}

double Particle::sphere_dUdr() const {
    if (norm(coord) < simulator::SPHERE_SIZE) return 0.0;
    return 2.0*(norm(coord) - simulator::SPHERE_SIZE);
    //return step_calculator::SPHERE_COEF;
}
} // smd

#endif