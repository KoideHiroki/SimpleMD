#ifndef SOAP_HPP
#define SOAP_HPP

#include "./constants.hpp"
#include "particle.hpp"
#include <random>

namespace smd {
class Soap {
public:
    Soap(const double init_x, const double init_y, const double init_z, std::mt19937& random_engine);
    Particle head;
    Particle tail;
    std::normal_distribution<double> head_dist;
    std::normal_distribution<double> tail_dist;
};

Soap::Soap(const double init_x, const double init_y, const double init_z, std::mt19937& random_engine)
    : head_dist(0.0, step_calculator::SOAP_HEAD_STD), tail_dist(0.0, step_calculator::SOAP_TAIL_STD), head(init_x-soap::SPRING_R0/2.0, init_y, init_z, soap::HEAD_WEIGHT, head_dist(random_engine), head_dist(random_engine), head_dist(random_engine)), tail(init_x+soap::SPRING_R0/2.0, init_y, init_z, soap::TAIL_WEIGHT, tail_dist(random_engine), tail_dist(random_engine), tail_dist(random_engine))
{}

} // smd

#endif