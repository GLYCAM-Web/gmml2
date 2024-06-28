#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <cmath>
#include <limits>
#include <iomanip> // setw

using cds::Coordinate;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
Coordinate::Coordinate(double x, double y, double z) : x_(x), y_(y), z_(z)
{}

Coordinate::Coordinate(const std::string x, const std::string y, const std::string z)
{
    try
    {
        x_ = std::stod(x);
        y_ = std::stod(y);
        z_ = std::stod(z);
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Could not convert these strings to doubles: " + x + ", " + y + ", " + z + ", ");
        throw;
    }
}

// Move Ctor.
Coordinate::Coordinate(Coordinate&& other) noexcept : Coordinate()
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Coordinate move ctor\n");
    x_ = std::move(other.x_);
    y_ = std::move(other.y_);
    z_ = std::move(other.z_);
}

Coordinate::Coordinate(const Coordinate& coordinate) : x_(coordinate.x_), y_(coordinate.y_), z_(coordinate.z_)
{
    // gmml::log(__LINE__,__FILE__,gmml::INF, "Coordinate copy ctor for " + std::to_string(x_) + " " +
    // std::to_string(y_) + " " + std::to_string(z_));
}

Coordinate& Coordinate::operator=(const Coordinate& other)
{
    x_ = other.x_;
    y_ = other.y_;
    z_ = other.z_;
    // gmml::log(__LINE__,__FILE__,gmml::INF, "Assign copy ctor for " + std::to_string(x_) + " " + std::to_string(y_) +
    // " " + std::to_string(z_));
    return *this;
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void Coordinate::Translate(double x, double y, double z)
{
    x_ += x;
    y_ += y;
    z_ += z;
}

double Coordinate::Distance(const Coordinate* coordinate) const
{
    double dist = (x_ - coordinate->GetX()) * (x_ - coordinate->GetX()) +
                  (y_ - coordinate->GetY()) * (y_ - coordinate->GetY()) +
                  (z_ - coordinate->GetZ()) * (z_ - coordinate->GetZ());
    if (dist > 0.00000001) // can sometimes measure distance to self, in which case get sqrt(0), which should be fine
                           // but zero is funky and somtimes is actually slightly negative.
    {
        return sqrt(dist);
    }
    return 0.0;
}

double Coordinate::length() const
{
    return sqrt((this->GetX() * this->GetX()) + (this->GetY() * this->GetY()) + (this->GetZ() * this->GetZ()));
}

void Coordinate::Normalize()
{
    double length = this->length();
    if (length != 0.0)
    {
        x_ = x_ / length;
        y_ = y_ / length;
        z_ = z_ / length;
    }
}

double Coordinate::DotProduct(const Coordinate& coordinate)
{
    return ((x_ * coordinate.x_) + (y_ * coordinate.y_) + (z_ * coordinate.z_));
}

// Should this not return a coord instead of altering this one?
void Coordinate::CrossProduct(const Coordinate& coordinate)
{
    double x = x_;
    double y = y_;
    double z = z_;
    x_       = (y * coordinate.z_) - (coordinate.y_ * z);
    y_       = (z * coordinate.x_) - (coordinate.z_ * x);
    z_       = (x * coordinate.y_) - (coordinate.x_ * y);
}

// These can all be better, the call is weird:
void Coordinate::operator+(const Coordinate& coordinate)
{
    x_ += coordinate.x_;
    y_ += coordinate.y_;
    z_ += coordinate.z_;
}

void Coordinate::operator+(const double addition)
{
    x_ += addition;
    y_ += addition;
    z_ += addition;
}

void Coordinate::operator-(const Coordinate& coordinate)
{
    x_ -= coordinate.x_;
    y_ -= coordinate.y_;
    z_ -= coordinate.z_;
}

void Coordinate::operator/(const Coordinate& coordinate)
{
    x_ /= coordinate.x_;
    y_ /= coordinate.y_;
    z_ /= coordinate.z_;
}

void Coordinate::operator/(const double divisor)
{
    x_ /= divisor;
    y_ /= divisor;
    z_ /= divisor;
}

void Coordinate::operator*(const double multiplier)
{
    x_ *= multiplier;
    y_ *= multiplier;
    z_ *= multiplier;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void Coordinate::Print(std::ostream& out) const
{
    if (this->GetX() == constants::dNotSet || this->GetY() == constants::dNotSet || this->GetZ() == constants::dNotSet)
    {
        out << std::setw(10) << " "
            << ", " << std::setw(10) << " "
            << ", " << std::setw(10) << " ";
    }
    else
    {
        out << std::setw(10) << x_ << ", " << std::setw(10) << y_ << ", " << std::setw(10) << z_;
    }
}

std::string Coordinate::ToString() const
{
    std::stringstream ss;
    this->Print(ss);
    return ss.str();
}
