#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"

#include <cmath>
#include <array>
#include <vector>

namespace
{
    using cds::Coordinate;
    using cds::Sphere;

    double radius(const Sphere& sphere)
    {
        return sphere.radius;
    }

    Coordinate center(const Sphere& sphere)
    {
        return sphere.center;
    }

    Sphere boundingSphereInitialEstimate(const std::vector<Sphere>& spheres)
    {
        const Coordinate init = center(spheres[0]);
        double d              = radius(spheres[0]);
        double x              = init.nth(0);
        double y              = init.nth(1);
        double z              = init.nth(2);
        std::array<double, 3> minValue {x - d, y - d, z - d};
        std::array<double, 3> maxValue {x + d, y + d, z + d};
        std::array<size_t, 3> minId {0, 0, 0};
        std::array<size_t, 3> maxId {0, 0, 0};

        for (size_t k = 1; k < spheres.size(); k++)
        {
            auto& a       = spheres[k];
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
            auto& minp  = spheres[minId[n]];
            auto& maxp  = spheres[maxId[n]];
            double span = cds::distance(center(minp), center(maxp)) + radius(minp) + radius(maxp);
            if (span > maxSpan)
            {
                maxSpan  = span;
                maxSpanN = n;
            }
        }
        auto finalCenter = cds::scaleBy(0.5, center(spheres[minId[maxSpanN]]) + center(spheres[maxId[maxSpanN]]));
        return Sphere {0.5 * maxSpan, finalCenter};
    }

    Sphere boundingSphereIncludingAllSpheres(Sphere sphere, const std::vector<Sphere>& spheres)
    {
        for (auto& a : spheres)
        {
            sphere = cds::boundingSphereIncluding(sphere, a);
        }
        return sphere;
    }
} // namespace

cds::Sphere cds::boundingSphereIncluding(Sphere sphere, const Sphere include)
{
    Coordinate diff   = center(include) - sphere.center;
    double diffLength = length(diff);
    double dist       = diffLength + radius(include);
    if (dist > sphere.radius)
    {
        double newRadius     = 0.5 * (sphere.radius + dist);
        Coordinate newCenter = sphere.center + cds::scaleBy((dist - sphere.radius) / (2.0 * diffLength), diff);
        sphere               = Sphere {newRadius, newCenter};
    }
    return sphere;
}

cds::Sphere cds::boundingSphere(const std::vector<Sphere>& spheres)
{
    Sphere sphere = boundingSphereInitialEstimate(spheres);
    return boundingSphereIncludingAllSpheres(sphere, spheres);
}
