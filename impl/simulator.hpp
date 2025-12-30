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
    Simulator();
    void run();
private:
    void step();
    void write_log(std::ofstream& out, const int loop_idx);
    void init_waters();
    void init_soaps();
    std::vector<Water> waters;
    std::vector<Soap> soaps;
    StepCalculator step_calculator;

    std::mt19937 random_engine;
    std::uniform_real_distribution<double> init_coord_dist;
};

Simulator::Simulator()
    : random_engine(simulator::SEED), init_coord_dist(-(simulator::SPHERE_SIZE-10.0), simulator::SPHERE_SIZE-10.0)
{
    init_waters();
    init_soaps();
}

void Simulator::run() {
    std::ofstream out(simulator::OUT_PATH);
    if (!out) {
        std::cerr << "cannot create: " << simulator::OUT_PATH << std::endl;
        std::exit(1);
    }
    out << "<meta>" << std::endl;
    out << "water_num " << simulator::WATER_NUM << std::endl;
    out << "soap_num " << simulator::SOAP_NUM << std::endl;
    out << "</meta>" << std::endl;
    for (int relax_idx = 0; relax_idx < simulator::RELAX_STEP_NUM; ++relax_idx) {
        const auto [relaxed_waters, relaxed_soaps] = step_calculator.relax(waters, soaps);
        waters = relaxed_waters;
        soaps = relaxed_soaps;
    }
    for (int loop_idx = 0; loop_idx < simulator::LOOP_NUM; ++loop_idx) {
        if (loop_idx % simulator::SAVE_STEP_NUM == 0) {
            write_log(out, loop_idx);
        }
        step();
    }
}

void Simulator::step() {
    const auto [new_waters, new_soaps] = step_calculator.calc(waters, soaps, random_engine);
    waters = new_waters;
    soaps = new_soaps;
}

void Simulator::write_log(std::ofstream& out, const int loop_idx) {
    out << "<frame " << loop_idx << ">" << std::endl;
    out << "<waters>" << std::endl;
    for (const auto& water: waters) {
        out << water.atom.coord[0] << " " << water.atom.coord[1] << " " << water.atom.coord[2] << std::endl;
    }
    out << "</waters>" << std::endl;
    out << "<soaps>" << std::endl;
    for (const auto& soap: soaps) {
        out << soap.head.coord[0] << " " << soap.head.coord[1] << " " << soap.head.coord[2] << " " << soap.tail.coord[0] << " " << soap.tail.coord[1] << " " << soap.tail.coord[2] << std::endl;
    }
    out << "</soaps>" << std::endl;
    out << "</frame>" << std::endl;
    std::cout << "step: " << loop_idx << std::endl;
}

void Simulator::init_waters() {
    for (int water_idx = 0; water_idx < simulator::WATER_NUM; ++water_idx) {
        waters.push_back(Water(init_coord_dist(random_engine), init_coord_dist(random_engine), init_coord_dist(random_engine), random_engine));
    }
}

void Simulator::init_soaps() {
    for (int soap_idx = 0; soap_idx < simulator::SOAP_NUM; ++soap_idx) {
        soaps.push_back(Soap(init_coord_dist(random_engine), init_coord_dist(random_engine), init_coord_dist(random_engine), random_engine));
    }
}

} // smd

#endif