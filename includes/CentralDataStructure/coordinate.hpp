#ifndef INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP

#include <iostream>
#include <vector>

namespace cds
{
    class Coordinate
    {
      public:
        Coordinate(double x, double y, double z) : x_(x), y_(y), z_(z)
        {}

        Coordinate(const Coordinate& a) : Coordinate(a.GetX(), a.GetY(), a.GetZ()) {};
        Coordinate& operator=(const Coordinate&);

        //////////////////////////////////////////////////////////
        //                           ACCESSOR                   //
        //////////////////////////////////////////////////////////
        double GetX() const
        {
            return x_;
        }

        double GetY() const
        {
            return y_;
        }

        double GetZ() const
        {
            return z_;
        }

        //////////////////////////////////////////////////////////
        //                           MUTATOR                    //
        //////////////////////////////////////////////////////////
        void SetX(const double x)
        {
            x_ = x;
        }

        void SetY(const double y)
        {
            y_ = y;
        }

        void SetZ(const double z)
        {
            z_ = z;
        }

        void Set(const double x, const double y, const double z)
        {
            x_ = x;
            y_ = y;
            z_ = z;
        }

        //////////////////////////////////////////////////////////
        //                         FUNCTIONS                    //
        //////////////////////////////////////////////////////////
        void Translate(const double x, const double y, const double z);
        double Distance(const Coordinate* coordinate) const;
        double length() const;
        void Normalize();
        double DotProduct(const Coordinate& coordinate);
        void CrossProduct(const Coordinate& coordinate);
        void operator+(const Coordinate& coordinate);
        void operator+(const double addition);
        void operator-(const Coordinate& coordinate);
        void operator/(const Coordinate& coordinate);
        void operator/(const double divisor);
        void operator*(const double multiplier);

        bool operator==(const Coordinate& rhs) const
        {
            return (this->GetX() == rhs.GetX() && this->GetY() == rhs.GetY() && this->GetZ() == rhs.GetZ());
        }

        Coordinate& operator+=(const Coordinate& rhs)
        {
            x_ += rhs.GetX();
            y_ += rhs.GetY();
            z_ += rhs.GetZ();
            return *this;
        }

        // Coordinate operator-(const Coordinate& rhs) const { return Coordinate(x_ - rhs.x_, y_ - rhs.y_, z_ -
        // rhs.z_);}
        //////////////////////////////////////////////////////////
        //                     DISPLAY FUNCTIONS                //
        //////////////////////////////////////////////////////////
        void Print(std::ostream& out = std::cerr) const;
        std::string ToString() const;

      private:
        //////////////////////////////////////////////////////////
        //                         ATTRIBUTES                   //
        //////////////////////////////////////////////////////////
        double x_;
        double y_;
        double z_;
    };

    inline bool withinDistance(double distance, const Coordinate& a, const Coordinate& b)
    {
        double dx = a.GetX() - b.GetX();
        double dy = a.GetY() - b.GetY();
        double dz = a.GetZ() - b.GetZ();
        return (dx * dx + dy * dy + dz * dz) < distance * distance;
    }

    Coordinate coordinateFromStrings(const std::string& x, const std::string& y, const std::string& z);
} // namespace cds
#endif
