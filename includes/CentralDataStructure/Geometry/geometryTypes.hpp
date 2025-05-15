#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_GEOMETRYTYPES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_GEOMETRYTYPES_HPP

#include <array>
#include <vector>

namespace cds
{
    class Coordinate
    {
      public:
        Coordinate(double x, double y, double z) : values_({x, y, z})
        {}

        inline double nth(int n) const
        {
            return values_[n];
        }

        inline double GetX() const
        {
            return nth(0);
        }

        inline double GetY() const
        {
            return nth(1);
        }

        inline double GetZ() const
        {
            return nth(2);
        }

        inline Coordinate operator+(const Coordinate& a) const
        {
            auto add = [&](int n)
            {
                return nth(n) + a.nth(n);
            };
            return {add(0), add(1), add(2)};
        }

        inline Coordinate operator-(const Coordinate& a) const
        {
            auto sub = [&](int n)
            {
                return nth(n) - a.nth(n);
            };
            return {sub(0), sub(1), sub(2)};
        }

      private:
        std::array<double, 3> values_;
    };

    struct Sphere
    {
        double radius;
        Coordinate center;
    };
} // namespace cds

#endif
