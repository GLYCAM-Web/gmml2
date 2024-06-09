#ifndef INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP

#include <array>
#include <iostream>

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

        void Print(std::ostream& out = std::cerr) const;
        std::string ToString() const;

      private:
        std::array<double, 3> values_;
    };

    inline double squaredLength(const Coordinate& a)
    {
        auto sq = [&](int n)
        {
            double d = a.nth(n);
            return d * d;
        };
        return sq(0) + sq(1) + sq(2);
    }

    inline double squaredDistance(const Coordinate& a, const Coordinate& b)
    {
        return squaredLength(a - b);
    }

    inline double dotProduct(const Coordinate& a, const Coordinate& b)
    {
        auto dot = [&](int n)
        {
            return a.nth(n) * b.nth(n);
        };
        return dot(0) + dot(1) + dot(2);
    }

    inline bool withinDistance(double distance, const Coordinate& a, const Coordinate& b)
    {
        return squaredDistance(a, b) < distance * distance;
    }

    double length(const Coordinate& a);
    double distance(const Coordinate& a, const Coordinate& b);
    Coordinate scaleBy(double factor, const Coordinate& a);
    Coordinate normal(const Coordinate& a);
    Coordinate crossProduct(const Coordinate& a, const Coordinate& b);

    Coordinate coordinateFromStrings(const std::string& x, const std::string& y, const std::string& z);
} // namespace cds
#endif
