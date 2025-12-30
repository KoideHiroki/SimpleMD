#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <string>

namespace smd {
    namespace step_calculator {
        const double SOFT_REPULSIVE_D = 3.0;
        const double RELAX_COEF = 1.0;
        const double RELAX_SPRING_COEF = 0.1;
        const double GAMMA = 1.0;
        const double KBT = 3.0;
        const double DT = 0.1;
        const double EXCLUDED_D = 0.9;

        const double WATER_EPSILON = 2.0;
        const double WATER_SIGMA = 0.9;
        const double WATER_HEAD_EPSILON = 1.2;
        const double WATER_HEAD_SIGMA = 0.9;
        const double HEAD_HEAD_EPSILON = 0.3;
        const double HEAD_HEAD_SIGMA = 0.9;
        const double TAIL_TAIL_EPSILON = 0.8;
        const double TAIL_TAIL_SIGMA = 0.9;
        const double HEAD_TAIL_COEF = 3.0;
        const double WATER_TAIL_COEF = 1.0;

        const double SPHERE_COEF = 1.0;
    }

    namespace particle {
        const double SOFT_REPULSIVE_A = 0.5;
        const double REPULSIVE_D = 3.0;
        const double FMAX = 50.0;
    }

    namespace soap {
        const double SPRING_K = 1.5;
        const double SPRING_R0 = 2.0;
        const double HEAD_WEIGHT = 1.0;
        const double TAIL_WEIGHT = 1.0;
    }

    namespace water {
        const double WEIGHT = 1.0;
    }

    namespace simulator {
        const double SPHERE_SIZE = 50.0;
        const int RELAX_STEP_NUM = 1000;
        const int WATER_NUM = 300;
        const int SOAP_NUM = 200;
        const uint32_t SEED = 12345678;
        const int LOOP_NUM = 100001;
        const int SAVE_STEP_NUM = 100;
        const std::string OUT_PATH = "exe.log";

    }

    namespace step_calculator {
        const double WATER_STD = ((2.0*GAMMA*KBT)/water::WEIGHT)*DT;
        const double SOAP_HEAD_STD = ((2.0*GAMMA*KBT)/soap::HEAD_WEIGHT)*DT;
        const double SOAP_TAIL_STD = ((2.0*GAMMA*KBT)/soap::TAIL_WEIGHT)*DT;
    }
} // smd


#endif