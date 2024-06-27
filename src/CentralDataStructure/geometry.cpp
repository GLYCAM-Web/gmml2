#include "includes/CentralDataStructure/geometry.hpp"

#include "includes/CentralDataStructure/coordinate.hpp"

#include <cmath>
#include <array>
#include <vector>

namespace
{
    using cds::Coordinate;

    cds::Sphere boundingSphereInitialEstimate(const std::vector<Coordinate>& points)
    {
        const Coordinate& init = points[0];
        double x               = init.nth(0);
        double y               = init.nth(1);
        double z               = init.nth(2);
        std::array<double, 3> minValue {x, y, z};
        std::array<double, 3> maxValue {x, y, z};
        std::array<size_t, 3> minId {0, 0, 0};
        std::array<size_t, 3> maxId {0, 0, 0};

        for (size_t k = 1; k < points.size(); k++)
        {
            auto& a = points[k];
            for (size_t n = 0; n < 3; n++)
            {
                double nth = a.nth(n);
                if (nth < minValue[n])
                {
                    minValue[n] = nth;
                    minId[n]    = k;
                }
                if (nth > maxValue[n])
                {
                    maxValue[n] = nth;
                    maxId[n]    = k;
                }
            }
        }
        size_t maxSpanN = 0;
        double maxSpan  = 0.0;
        for (size_t n = 0; n < 3; n++)
        {
            double span = cds::squaredDistance(points[minId[n]], points[maxId[n]]);
            if (span > maxSpan)
            {
                maxSpan  = span;
                maxSpanN = n;
            }
        }

        return cds::Sphere {0.5 * std::sqrt(maxSpan),
                            cds::scaleBy(0.5, points[minId[maxSpanN]] + points[maxId[maxSpanN]])};
    }

    cds::Sphere boundingSphereIncludingAllPoints(cds::Sphere sphere, const std::vector<Coordinate>& points)
    {
        for (auto& a : points)
        {
            Coordinate diff = a - sphere.center;
            double sqDist   = squaredLength(diff);
            if (sqDist > (sphere.radius * sphere.radius))
            {
                double dist          = std::sqrt(sqDist);
                double newRadius     = 0.5 * (sphere.radius + dist);
                Coordinate newCenter = sphere.center + cds::scaleBy((dist - sphere.radius) / (2.0 * dist), diff);
                sphere               = cds::Sphere {newRadius, newCenter};
            }
        }
        return sphere;
    }
} // namespace

cds::Sphere cds::boundingSphere(const std::vector<Coordinate>& points)
{
    if (points.empty())
    {
        return Sphere {
            0.0, {0.0, 0.0, 0.0}
        };
    }
    else
    {
        Sphere sphere = boundingSphereInitialEstimate(points);
        return boundingSphereIncludingAllPoints(sphere, points);
    }
}
