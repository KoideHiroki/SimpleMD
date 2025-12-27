#ifndef SOAP_HPP
#define SOAP_HPP

#include "./constants.hpp"

namespace smd {
class Soap {
public:
    Soap(const double init_x, const double init_y);
    double head_x;
    double head_y;
    double tail_x;
    double tail_y;
};

Soap::Soap(const double init_x, const double init_y) {
    head_x = init_x-soap::MOL_LEN/2.0;
    head_y = init_y;
    tail_x = init_x+soap::MOL_LEN/2.0;
    tail_y = init_y;
}
} // smd

#endif