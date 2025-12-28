#ifndef STEP_CALCULATOR_HPP
#define STEP_CALCULATOR_HPP

#include "./water.hpp"
#include "./soap.hpp"
#include <tuple>
#include <vector>
#include <random>

namespace smd {
class StepCalculator {
public:
    StepCalculator() {}
    std::tuple<std::vector<Water>,std::vector<Soap>> calc(const std::vector<Water>& waters, const std::vector<Soap>& soaps, std::mt19937& random_engine);
private:
};

std::tuple<std::vector<Water>,std::vector<Soap>> StepCalculator::calc(const std::vector<Water>& waters, const std::vector<Soap>& soaps, std::mt19937& random_engine) {
    return std::make_tuple(waters, soaps);
}
} // smd

#endif