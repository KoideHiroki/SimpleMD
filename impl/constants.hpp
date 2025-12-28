#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace smd {
    namespace step_calculator {
        const double BOX_SIZE = 100.0;
        const double REPULSIVE_D = 3.0;
        const double RELAX_COEF = 0.1;
        const double RELAX_SPRING_COEF = 0.1;
        const double GAMMA = 1.0;
        const double KBT = 1.0;
        const double DT = 0.01;

        const double WATER_EPSILON = 1.0;
        const double WATER_SIGMA = 1.0;
        const double WATER_HEAD_EPSILON = 1.0;
        const double WATER_HEAD_SIGMA = 1.0;
        const double HEAD_HEAD_EPSILON = 1.0;
        const double HEAD_HEAD_SIGMA = 1.0;
        const double TAIL_TAIL_EPSILON = 1.0;
        const double TAIL_TAIL_SIGMA = 1.0;
        const double HEAD_TAIL_COEF = 1.0;
        const double WATER_TAIL_COEF = 1.0;
    }

    namespace particle {
        const double SOFT_REPULSIVE_A = 0.1;
    }

    namespace soap {
        const double SPRING_K = 1.0;
        const double SPRING_R0 = 2.0;
        const double HEAD_WEIGHT = 1.0;
        const double TAIL_WEIGHT = 1.0;
    }

    namespace water {
        const double WEIGHT = 1.0;
    }
} // smd


#endif