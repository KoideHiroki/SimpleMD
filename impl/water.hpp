#ifndef WATER_HPP
#define WATER_HPP

#include "particle.hpp"

namespace smd {
class Water {
public:
    Water(const double init_x, const double init_y, const double init_z);
    Particle atom;
};

Water::Water(const double init_x, const double init_y, const double init_z)
    : atom(init_x, init_y, init_z)
{}
} // smd

#endif