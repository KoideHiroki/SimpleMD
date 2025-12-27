#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include <iostream>
#include <fstream>
#include <string>

namespace smd {
class Simulator {
public:
    Simulator();
    void run(const int loop_num, const int save_step_num, const std::string& out_path);
private:
    void step();
    void write_log(std::ofstream& out)
};

Simulator::Simulator() {

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

}

void Simulator::write_log(std::ofstream& out) {

}
} // smd

#endif