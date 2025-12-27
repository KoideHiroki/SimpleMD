#include "./impl/simulator.hpp"

int main() {
    smd::Simulator simulator(12345678);
    simulator.run(10, 1, "test.log");
}