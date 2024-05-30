#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <cmath>
#include <iomanip> // setw

using cds::Coordinate;

double cds::length(const Coordinate& a)
{
    return std::sqrt(squaredLength(a));
}

double cds::distance(const Coordinate& a, const Coordinate& b)
{
    return std::sqrt(squaredDistance(a, b));
}

Coordinate cds::scaleBy(double factor, const Coordinate& a)
{
    auto scaled = [&](int n)
    {
        return factor * a.nth(n);
    };
    return {scaled(0), scaled(1), scaled(2)};
}

Coordinate cds::normal(const Coordinate& a)
{
    return scaleBy(1.0 / length(a), a);
}

Coordinate cds::crossProduct(const Coordinate& a, const Coordinate& b)
{
    auto cross = [&](int n, int k)
    {
        return a.nth(n) * b.nth(k) - a.nth(k) * b.nth(n);
    };
    return {cross(1, 2), cross(2, 0), cross(0, 1)};
}

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
        out << std::setw(10) << GetX() << ", " << std::setw(10) << GetY() << ", " << std::setw(10) << GetZ();
    }
}

std::string Coordinate::ToString() const
{
    std::stringstream ss;
    this->Print(ss);
    return ss.str();
}

Coordinate cds::coordinateFromStrings(const std::string& x, const std::string& y, const std::string& z)
{
    try
    {
        return {std::stod(x), std::stod(y), std::stod(z)};
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Could not convert these strings to doubles: " + x + ", " + y + ", " + z + ", ");
        throw;
    }
}
