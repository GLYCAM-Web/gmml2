#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"

#include <cmath>
#include <array>
#include <vector>

namespace
{
    using cds::Coordinate;

    double radius(const Coordinate&)
    {
        return 0.0;
    }

    double radius(const cds::Sphere& sphere)
    {
        return sphere.radius;
    }

    Coordinate center(const Coordinate& coord)
    {
        return coord;
    }

    Coordinate center(const cds::Sphere& sphere)
    {
        return sphere.center;
    }

    template<class T> cds::Sphere boundingSphereInitialEstimate(const std::vector<T>& points)
    {
        const Coordinate init = center(points[0]);
        double d              = radius(points[0]);
        double x              = init.nth(0);
        double y              = init.nth(1);
        double z              = init.nth(2);
        std::array<double, 3> minValue {x - d, y - d, z - d};
        std::array<double, 3> maxValue {x + d, y + d, z + d};
        std::array<size_t, 3> minId {0, 0, 0};
        std::array<size_t, 3> maxId {0, 0, 0};

        for (size_t k = 1; k < points.size(); k++)
        {
            auto& a       = points[k];
            Coordinate pt = center(a);
            double r      = radius(a);
            for (size_t n = 0; n < 3; n++)
            {
                double nth = pt.nth(n);
                if (nth - r < minValue[n])
                {
                    minValue[n] = nth - r;
                    minId[n]    = k;
                }
                if (nth + r > maxValue[n])
                {
                    maxValue[n] = nth + r;
                    maxId[n]    = k;
                }
            }
        }
        size_t maxSpanN = 0;
        double maxSpan  = 0.0;
        for (size_t n = 0; n < 3; n++)
        {
            auto& minp  = points[minId[n]];
            auto& maxp  = points[maxId[n]];
            double span = cds::distance(center(minp), center(maxp)) + radius(minp) + radius(maxp);
            if (span > maxSpan)
            {
                maxSpan  = span;
                maxSpanN = n;
            }
        }

        return cds::Sphere {0.5 * maxSpan,
                            cds::scaleBy(0.5, center(points[minId[maxSpanN]]) + center(points[maxId[maxSpanN]]))};
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

    cds::Sphere boundingSphereIncludingAllSpheres(cds::Sphere sphere, const std::vector<cds::Sphere>& spheres)
    {
        for (auto& a : spheres)
        {
            Coordinate diff   = center(a) - sphere.center;
            double diffLength = length(diff);
            double dist       = diffLength + radius(a);
            if (dist > sphere.radius)
            {
                double newRadius     = 0.5 * (sphere.radius + dist);
                Coordinate newCenter = sphere.center + cds::scaleBy((dist - sphere.radius) / (2.0 * diffLength), diff);
                sphere               = cds::Sphere {newRadius, newCenter};
            }
        }
        return sphere;
    }
} // namespace

cds::Sphere cds::boundingSphere(const std::vector<Coordinate>& points)
{
    Sphere sphere = boundingSphereInitialEstimate(points);
    return boundingSphereIncludingAllPoints(sphere, points);
}

cds::Sphere cds::boundingSphere(const std::vector<Sphere>& spheres)
{
    Sphere sphere = boundingSphereInitialEstimate(spheres);
    return boundingSphereIncludingAllSpheres(sphere, spheres);
}
