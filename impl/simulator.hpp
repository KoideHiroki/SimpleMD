#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "constants.hpp"
#include "./water.hpp"
#include "./soap.hpp"
#include "./step_calculator.hpp"

namespace smd {
class Simulator {
public:
    Simulator(const std::uint32_t seed);
    void run(const int loop_num, const int save_step_num, const std::string& out_path);
private:
    void step();
    void write_log(std::ofstream& out);
    void init_waters();
    void init_soaps();
    std::vector<Water> waters;
    std::vector<Soap> soaps;
    StepCalculator step_calculator;

    std::mt19937 random_engine;
    std::uniform_real_distribution<double> init_coord_dist;
};

Simulator::Simulator(const std::uint32_t seed)
    : random_engine(seed), init_coord_dist(-step_calculator::BOX_SIZE-10.0, step_calculator::BOX_SIZE-10.0)
{
    init_waters();
    init_soaps();
}

void Simulator::run(const int loop_num, const int save_step_num, const std::string& out_path) {
    std::ofstream out(out_path);
    if (!out) {
        std::cerr << "cannot create: " << out_path << std::endl;
        std::exit(1);
    }
    for (int loop_idx = 0; loop_idx < loop_num; ++loop_idx) {
        step();
        if (loop_idx % save_step_num == 0) {
            write_log(out);
        }
    }
}

void Simulator::step() {
    const auto [new_waters, new_soaps] = step_calculator.calc(waters, soaps, random_engine);
    waters = new_waters;
    soaps = new_soaps;
}

void Simulator::write_log(std::ofstream& out) {
}

void Simulator::init_waters() {
    for (int water_idx = 0; water_idx < 100; ++water_idx) {
        waters.push_back(Water(init_coord_dist(random_engine), init_coord_dist(random_engine), init_coord_dist(random_engine)));
    }
}

void Simulator::init_soaps() {
    for (int soap_idx = 0; soap_idx < 100; ++soap_idx) {
        soaps.push_back(Soap(init_coord_dist(random_engine), init_coord_dist(random_engine), init_coord_dist(random_engine)));
    }
}

} // smd

#endif