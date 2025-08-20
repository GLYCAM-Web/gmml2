#ifndef INCLUDE_GEOMETRY_GEOMETRYTYPES_HPP
#define INCLUDE_GEOMETRY_GEOMETRYTYPES_HPP

#include <array>
#include <vector>

namespace gmml
{
    class Coordinate
    {
      public:
        Coordinate(double x, double y, double z) : values({x, y, z}) {}

        inline double nth(int n) const { return values[n]; }

        inline Coordinate operator+(const Coordinate& a) const
        {
            auto add = [&](int n) { return nth(n) + a.nth(n); };
            return {add(0), add(1), add(2)};
        }

        inline Coordinate operator-(const Coordinate& a) const
        {
            auto sub = [&](int n) { return nth(n) - a.nth(n); };
            return {sub(0), sub(1), sub(2)};
        }

      private:
        std::array<double, 3> values;
    };

    struct Sphere
    {
        double radius;
        Coordinate center;
    };
} // namespace gmml

#endif
