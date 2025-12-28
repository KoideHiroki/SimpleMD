#ifndef SOAP_HPP
#define SOAP_HPP

#include "./constants.hpp"
#include "particle.hpp"

namespace smd {
class Soap {
public:
    Soap(const double init_x, const double init_y, const double init_z);
    Particle head;
    Particle tail;
};

Soap::Soap(const double init_x, const double init_y, const double init_z)
    : head(init_x-soap::SPRING_R0/2.0, init_y, init_z, soap::HEAD_WEIGHT), tail(init_x+soap::SPRING_R0/2.0, init_y, init_z, soap::TAIL_WEIGHT)
{}

} // smd

#endif