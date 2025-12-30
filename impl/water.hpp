#ifndef WATER_HPP
#define WATER_HPP

#include "particle.hpp"
#include <random>

namespace smd {
class Water {
public:
    Water(const double init_x, const double init_y, const double init_z, std::mt19937& random_engine);
    Particle atom;
    std::normal_distribution<double> water_dist;

};

Water::Water(const double init_x, const double init_y, const double init_z, std::mt19937& random_engine)
    : water_dist(0.0, step_calculator::WATER_STD), atom(init_x, init_y, init_z, water::WEIGHT, water_dist(random_engine), water_dist(random_engine), water_dist(random_engine))
{}
} // smd

#endif