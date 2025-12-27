#ifndef WATER_HPP
#define WATER_HPP

namespace smd {
class Water {
public:
    Water(const double init_x, const double init_y);
    double x;
    double y;
};

Water::Water(const double init_x, const double init_y) {
    x = init_x;
    y = init_y;
}
} // smd

#endif