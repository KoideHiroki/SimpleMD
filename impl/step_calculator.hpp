#ifndef STEP_CALCULATOR_HPP
#define STEP_CALCULATOR_HPP

#include "./water.hpp"
#include "./soap.hpp"
#include <tuple>
#include <vector>

namespace smd {
class StepCalculator {
public:
    StepCalculator() {}
    std::tuple<std::vector<Water>,std::vector<Soap>> calc(const std::vector<Water>& waters, const std::vector<Soap>& soaps);
private:
};

std::tuple<std::vector<Water>,std::vector<Soap>> StepCalculator::calc(const std::vector<Water>& waters, const std::vector<Soap>& soaps) {
    return std::make_tuple(waters, soaps);
}
} // smd

#endif