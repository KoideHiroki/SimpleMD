#ifndef STEP_CALCULATOR_HPP
#define STEP_CALCULATOR_HPP

#include "./water.hpp"
#include "./soap.hpp"
#include "./particle.hpp"
#include "constants.hpp"
#include <tuple>
#include <vector>
#include <random>

namespace smd {
class StepCalculator {
public:
    StepCalculator();
    std::tuple<std::vector<Water>,std::vector<Soap>> calc(const std::vector<Water>& waters, const std::vector<Soap>& soaps, std::mt19937& random_engine);
    std::tuple<std::vector<Water>,std::vector<Soap>> relax(const std::vector<Water>& waters, const std::vector<Soap>& soaps);
private:
    std::array<double,3> calc_watar_force(const int water_idx, const std::vector<Water>& waters, const std::vector<Soap>& soaps);
    std::array<double,3> calc_soap_head_force(const int soap_idx, const std::vector<Water>& waters, const std::vector<Soap>& soaps);
    std::array<double,3> calc_soap_tail_force(const int soap_idx, const std::vector<Water>& waters, const std::vector<Soap>& soaps);

    std::normal_distribution<double> water_dist;
    std::normal_distribution<double> soap_head_dist;
    std::normal_distribution<double> soap_tail_dist;

    //const double water_std = ((2.0*step_calculator::GAMMA*step_calculator::KBT)/water::WEIGHT)/step_calculator::DT;
    //const double soap_head_std = ((2.0*step_calculator::GAMMA*step_calculator::KBT)/soap::HEAD_WEIGHT)/step_calculator::DT;
    //const double soap_tail_std = ((2.0*step_calculator::GAMMA*step_calculator::KBT)/soap::TAIL_WEIGHT)/step_calculator::DT;
};

StepCalculator::StepCalculator()
    : water_dist(0.0, step_calculator::WATER_STD), soap_head_dist(0.0, step_calculator::SOAP_HEAD_STD), soap_tail_dist(0.0, step_calculator::SOAP_TAIL_STD)
    //: water_dist(0.0, ((2.0*step_calculator::GAMMA*step_calculator::KBT)/water::WEIGHT)/step_calculator::DT),
    //    soap_head_dist(0.0, ((2.0*step_calculator::GAMMA*step_calculator::KBT)/soap::HEAD_WEIGHT)/step_calculator::DT),
    //    soap_tail_dist(0.0, ((2.0*step_calculator::GAMMA*step_calculator::KBT)/soap::TAIL_WEIGHT)/step_calculator::DT)
{}

std::tuple<std::vector<Water>,std::vector<Soap>> StepCalculator::calc(const std::vector<Water>& waters, const std::vector<Soap>& soaps, std::mt19937& random_engine) {
    const double coef = 1.0 - (step_calculator::GAMMA*step_calculator::DT)/2.0;
    auto new_waters = waters;
    auto new_soaps = soaps;

    std::vector<std::array<double,3>> water_forces_old;
    for (int water_idx = 0; water_idx < waters.size(); ++water_idx) {
        const auto f = calc_watar_force(water_idx, waters, soaps);
        water_forces_old.push_back(f);
    }

    std::vector<std::array<double,3>> soap_head_forces_old;
    for (int soap_idx = 0; soap_idx < soaps.size(); ++soap_idx) {
        const auto f = calc_soap_head_force(soap_idx, waters, soaps);
        soap_head_forces_old.push_back(f);
    }

    std::vector<std::array<double,3>> soap_tail_forces_old;
    for (int soap_idx = 0; soap_idx < soaps.size(); ++soap_idx) {
        const auto f = calc_soap_tail_force(soap_idx, waters, soaps);
        soap_tail_forces_old.push_back(f);
    }
    for (int water_idx = 0; water_idx < waters.size(); ++water_idx) {
        //std::cout << waters[water_idx].atom.random_memory[0] << std::endl;
        new_waters[water_idx].atom.coord += (step_calculator::DT*coef)*waters[water_idx].atom.velo + (0.5*step_calculator::DT*step_calculator::DT)*((1.0/water::WEIGHT)*water_forces_old[water_idx] + waters[water_idx].atom.random_memory);
    }

    for (int soap_idx = 0; soap_idx < soaps.size(); ++soap_idx) {
        new_soaps[soap_idx].head.coord += (step_calculator::DT*coef)*soaps[soap_idx].head.velo + (0.5*step_calculator::DT*step_calculator::DT)*((1.0/soap::HEAD_WEIGHT)*soap_head_forces_old[soap_idx] + soaps[soap_idx].head.random_memory);
        new_soaps[soap_idx].tail.coord += (step_calculator::DT*coef)*soaps[soap_idx].tail.velo + (0.5*step_calculator::DT*step_calculator::DT)*((1.0/soap::TAIL_WEIGHT)*soap_tail_forces_old[soap_idx] + soaps[soap_idx].tail.random_memory);
    }

    for (int water_idx = 0; water_idx < waters.size(); ++water_idx) {
        const std::array<double,3> new_random_memory = {water_dist(random_engine), water_dist(random_engine), water_dist(random_engine)};
        new_waters[water_idx].atom.velo = (coef * (coef + std::pow((step_calculator::GAMMA*step_calculator::DT)/2.0, 2.0)))*waters[water_idx].atom.velo + (step_calculator::DT/2.0)*((1.0/water::WEIGHT)*water_forces_old[water_idx]+(1.0/water::WEIGHT)*calc_watar_force(water_idx, new_waters, new_soaps) + waters[water_idx].atom.random_memory + new_random_memory);
        new_waters[water_idx].atom.random_memory = new_random_memory;
    }

    for (int soap_idx = 0; soap_idx < soaps.size(); ++soap_idx) {
        const std::array<double,3> new_head_random_memory = {soap_head_dist(random_engine), soap_head_dist(random_engine), soap_head_dist(random_engine)};
        const std::array<double,3> new_tail_random_memory = {soap_tail_dist(random_engine), soap_tail_dist(random_engine), soap_tail_dist(random_engine)};
        new_soaps[soap_idx].head.velo = (coef * (coef + std::pow((step_calculator::GAMMA*step_calculator::DT)/2.0, 2.0)))*soaps[soap_idx].head.velo + (step_calculator::DT/2.0)*((1.0/soap::HEAD_WEIGHT)*soap_head_forces_old[soap_idx]+(1.0/soap::HEAD_WEIGHT)*calc_soap_head_force(soap_idx, new_waters, new_soaps) + soaps[soap_idx].head.random_memory + new_head_random_memory);
        new_soaps[soap_idx].tail.velo = (coef * (coef + std::pow((step_calculator::GAMMA*step_calculator::DT)/2.0, 2.0)))*soaps[soap_idx].tail.velo + (step_calculator::DT/2.0)*((1.0/soap::TAIL_WEIGHT)*soap_tail_forces_old[soap_idx]+(1.0/soap::TAIL_WEIGHT)*calc_soap_tail_force(soap_idx, new_waters, new_soaps) + soaps[soap_idx].tail.random_memory + new_tail_random_memory);
        new_soaps[soap_idx].head.random_memory = new_head_random_memory;
        new_soaps[soap_idx].tail.random_memory = new_tail_random_memory;
    }
    //const auto diff = new_waters[0].atom.coord - waters[0].atom.coord;
    //std::cout << norm(diff) << std::endl;
    //const auto diff_n = norm(diff);
    //std::cout << diff_n << std::endl;
    //if (norm(diff) > 1.0) {
    //    std::cerr << "diverse position diff." << std::endl;
    //    std::exit(1);
    //}
    return std::make_tuple(new_waters, new_soaps);
}

std::tuple<std::vector<Water>,std::vector<Soap>> StepCalculator::relax(const std::vector<Water>& waters, const std::vector<Soap>& soaps) {
    auto new_waters = waters;
    auto new_soaps = soaps;
    for (int water_idx = 0; water_idx < waters.size(); ++water_idx) {
        std::array<double,3> force = {0.0, 0.0, 0.0};
        for (int other_water_idx = 0; other_water_idx < waters.size(); ++other_water_idx) {
            if (water_idx == other_water_idx) continue;
            const auto f = waters[water_idx].atom.calc_soft_repulsive(waters[other_water_idx].atom, step_calculator::SOFT_REPULSIVE_D);
            force += f;
        }
        for (const auto& soap : soaps) {
            const auto f_head = waters[water_idx].atom.calc_soft_repulsive(soap.head, step_calculator::SOFT_REPULSIVE_D);
            const auto f_tail = waters[water_idx].atom.calc_soft_repulsive(soap.tail, step_calculator::SOFT_REPULSIVE_D);
            force += f_head;
            force += f_tail;
        }
        //force = force / (norm(force)+1e-6);
        if (norm(waters[water_idx].atom.coord) > simulator::SPHERE_SIZE) {
            force += -step_calculator::SPHERE_COEF*(waters[water_idx].atom.coord/norm(waters[water_idx].atom.coord));
        }
        const auto diff = step_calculator::RELAX_COEF*force;
        new_waters[water_idx].atom.coord += diff;
    }

    for (int soap_idx = 0; soap_idx < soaps.size(); ++soap_idx) {
        std::array<double,3> force_head = {0.0, 0.0, 0.0};
        std::array<double,3> force_tail = {0.0, 0.0, 0.0};
        for (const auto& water : waters) {
            const auto f_head = soaps[soap_idx].head.calc_soft_repulsive(water.atom, step_calculator::SOFT_REPULSIVE_D);
            const auto f_tail = soaps[soap_idx].tail.calc_soft_repulsive(water.atom, step_calculator::SOFT_REPULSIVE_D);
            force_head += f_head;
            force_tail += f_tail;
        }
        for (int other_soap_idx = 0; other_soap_idx < soaps.size(); ++other_soap_idx) {
            if (soap_idx == other_soap_idx) continue;
            std::array<double,3> f_head = {0.0, 0.0, 0.0};
            std::array<double,3> f_tail = {0.0, 0.0, 0.0};
            f_head += soaps[soap_idx].head.calc_soft_repulsive(soaps[other_soap_idx].head, step_calculator::SOFT_REPULSIVE_D);
            f_head += soaps[soap_idx].head.calc_soft_repulsive(soaps[other_soap_idx].tail, step_calculator::SOFT_REPULSIVE_D);
            f_tail += soaps[soap_idx].tail.calc_soft_repulsive(soaps[other_soap_idx].head, step_calculator::SOFT_REPULSIVE_D);
            f_tail += soaps[soap_idx].tail.calc_soft_repulsive(soaps[other_soap_idx].tail, step_calculator::SOFT_REPULSIVE_D);
            force_head += f_head;
            force_tail += f_tail;
        }
        //force_head = force_head / (norm(force_head)+1e-6);
        //force_tail = force_tail / (norm(force_tail)+1e-6);
        if (norm(soaps[soap_idx].head.coord) > simulator::SPHERE_SIZE) {
            force_head += -step_calculator::SPHERE_COEF*(soaps[soap_idx].head.coord/norm(soaps[soap_idx].head.coord));
        }
        if (norm(soaps[soap_idx].tail.coord) > simulator::SPHERE_SIZE) {
            force_tail += -step_calculator::SPHERE_COEF*(soaps[soap_idx].tail.coord/norm(soaps[soap_idx].tail.coord));
        }
        const auto diff_head = step_calculator::RELAX_COEF*force_head;
        const auto diff_tail = step_calculator::RELAX_COEF*force_tail;
        new_soaps[soap_idx].head.coord += diff_head;
        new_soaps[soap_idx].tail.coord += diff_tail;

        const auto spring_force_head = soaps[soap_idx].head.calc_spring(soaps[soap_idx].tail, soap::SPRING_K, soap::SPRING_R0);
        const auto spring_force_tail = soaps[soap_idx].tail.calc_spring(soaps[soap_idx].head, soap::SPRING_K, soap::SPRING_R0);

        const auto spring_diff_head = step_calculator::RELAX_SPRING_COEF*spring_force_head;
        const auto spring_diff_tail = step_calculator::RELAX_SPRING_COEF*spring_force_tail;

        new_soaps[soap_idx].head.coord += spring_diff_head;
        new_soaps[soap_idx].tail.coord += spring_diff_tail;
    }
    return std::make_tuple(new_waters, new_soaps);
}

std::array<double,3> StepCalculator::calc_watar_force(const int water_idx, const std::vector<Water>& waters, const std::vector<Soap>& soaps) {
    std::array<double,3> force = {0.0, 0.0, 0.0};
    for (int other_water_idx = 0; other_water_idx < waters.size(); ++other_water_idx) {
        if (water_idx == other_water_idx) continue;
        force += waters[water_idx].atom.calc_LennardJones(waters[other_water_idx].atom, step_calculator::WATER_EPSILON, step_calculator::WATER_SIGMA);
        force += waters[water_idx].atom.calc_excluded(waters[other_water_idx].atom, step_calculator::EXCLUDED_D);
    }

    for (const auto& soap: soaps) {
        force += waters[water_idx].atom.calc_LennardJones(soap.head, step_calculator::WATER_HEAD_EPSILON, step_calculator::WATER_HEAD_SIGMA);
        force += waters[water_idx].atom.calc_repulsive(soap.tail, step_calculator::WATER_TAIL_COEF);
        //force += waters[water_idx].atom.calc_excluded(soap.tail, step_calculator::EXCLUDED_D);
        force += waters[water_idx].atom.calc_excluded(soap.head, step_calculator::EXCLUDED_D);
        force += waters[water_idx].atom.calc_excluded(soap.tail, step_calculator::EXCLUDED_D);
    }

    force += waters[water_idx].atom.calc_sphere();

    return force;
}

std::array<double,3> StepCalculator::calc_soap_head_force(const int soap_idx, const std::vector<Water>& waters, const std::vector<Soap>& soaps) {
    std::array<double,3> force = {0.0, 0.0, 0.0};
    for (const auto& water: waters) {
        force += soaps[soap_idx].head.calc_LennardJones(water.atom, step_calculator::WATER_HEAD_EPSILON, step_calculator::WATER_HEAD_SIGMA);
        force += soaps[soap_idx].head.calc_excluded(water.atom, step_calculator::EXCLUDED_D);
    }

    for (int other_soap_idx = 0; other_soap_idx < soaps.size(); ++other_soap_idx) {
        if (soap_idx == other_soap_idx) continue;
        force += soaps[soap_idx].head.calc_LennardJones(soaps[other_soap_idx].head, step_calculator::HEAD_HEAD_EPSILON, step_calculator::HEAD_HEAD_SIGMA);
        force += soaps[soap_idx].head.calc_repulsive(soaps[other_soap_idx].tail, step_calculator::HEAD_TAIL_COEF);
        //force += soaps[soap_idx].head.calc_excluded(soaps[other_soap_idx].tail, step_calculator::EXCLUDED_D);
        force += soaps[soap_idx].head.calc_excluded(soaps[other_soap_idx].head, step_calculator::EXCLUDED_D);
        force += soaps[soap_idx].head.calc_excluded(soaps[other_soap_idx].tail, step_calculator::EXCLUDED_D);
    }

    force += soaps[soap_idx].head.calc_spring(soaps[soap_idx].tail, soap::SPRING_K, soap::SPRING_R0);
    force += soaps[soap_idx].head.calc_sphere();

    return force;
}

std::array<double,3> StepCalculator::calc_soap_tail_force(const int soap_idx, const std::vector<Water>& waters, const std::vector<Soap>& soaps) {
    std::array<double,3> force = {0.0, 0.0, 0.0};
    for (const auto& water: waters) {
        force += soaps[soap_idx].tail.calc_repulsive(water.atom, step_calculator::WATER_TAIL_COEF);
        //force += soaps[soap_idx].tail.calc_excluded(water.atom, step_calculator::EXCLUDED_D);
        force += soaps[soap_idx].tail.calc_excluded(water.atom, step_calculator::EXCLUDED_D);
    }

    for (int other_soap_idx = 0; other_soap_idx < soaps.size(); ++other_soap_idx) {
        if (soap_idx == other_soap_idx) continue;
        force += soaps[soap_idx].tail.calc_repulsive(soaps[other_soap_idx].head, step_calculator::HEAD_TAIL_COEF);
        //force += soaps[soap_idx].tail.calc_excluded(soaps[other_soap_idx].head, step_calculator::EXCLUDED_D);
        force += soaps[soap_idx].tail.calc_LennardJones(soaps[other_soap_idx].tail, step_calculator::TAIL_TAIL_EPSILON, step_calculator::TAIL_TAIL_SIGMA);
        force += soaps[soap_idx].tail.calc_excluded(soaps[other_soap_idx].head, step_calculator::EXCLUDED_D);
        force += soaps[soap_idx].tail.calc_excluded(soaps[other_soap_idx].tail, step_calculator::EXCLUDED_D);
    }

    force += soaps[soap_idx].tail.calc_spring(soaps[soap_idx].head, soap::SPRING_K, soap::SPRING_R0);
    force += soaps[soap_idx].tail.calc_sphere();

    return force;
}
} // smd

#endif