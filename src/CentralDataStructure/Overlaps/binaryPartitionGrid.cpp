#include "includes/CentralDataStructure/Overlaps/binaryPartitionGrid.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"

#include <array>
#include <vector>

namespace
{
    size_t gridWidth(size_t subdivisions)
    {
        size_t width = 1;
        for (size_t n = 0; n < subdivisions; n++)
        {
            width *= 2;
        }
        return width;
    }

    size_t gridSize(size_t subdivisions)
    {
        size_t width = gridWidth(subdivisions);
        return width * width * width;
    }

    std::vector<int> gridNeighborOffsets(size_t subdivisions, const cds::GridStride& stride)
    {
        switch (subdivisions)
        {
            // loop will create duplicate entries if grid is smaller than 3x3x3, which is the case for 0-1
            case 0:
                return {0};
            case 1:
                return {-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
            default:
                std::array<int, 3> range {-1, 0, 1};
                std::vector<int> result;
                result.reserve(27);
                for (int z : range)
                {
                    for (int y : range)
                    {
                        for (int x : range)
                        {
                            result.push_back(x * stride[0] + y * stride[1] + z * stride[2]);
                        }
                    }
                }
                return result;
        }
    }
} // namespace

cds::BinaryPartitionGridBounds cds::binaryPartitionGridBounds(size_t subdivisions, double distance,
                                                              const Coordinate& coord)
{
    GridBoundsAxes limits;
    size_t width     = gridWidth(subdivisions);
    size_t halfWidth = width / 2;
    for (size_t n = 0; n < 3; n++)
    {
        limits[n].values.resize(width - 1);
        double c      = coord.nth(n);
        auto& values  = limits[n].values;
        int stride    = halfWidth;
        int m         = 1;
        size_t offset = 0;
        for (size_t sub = 0; sub < subdivisions; sub++)
        {
            for (int k = 0; k < m; k++)
            {
                int units          = (2 * k + 1) * stride - halfWidth;
                values[k + offset] = c + units * distance;
            }
            offset += m;
            stride /= 2;
            m      *= 2;
        }
    }
    GridStride strides {1, width, width * width};
    return BinaryPartitionGridBounds {subdivisions, strides, limits};
}

void cds::BinaryPartitionGrid::init(const BinaryPartitionGridBounds bounds, const std::vector<Coordinate>& coordinates)
{
    this->bounds_     = bounds;
    auto subdivisions = bounds.subdivisions;
    size_t size       = gridSize(subdivisions);
    auto& result      = this->grid_;
    result.resize(size);
    // special case to keep subdivision 0 relatively low-cost
    if (size == 1)
    {
        auto& target = result[0];
        target.clear();
        target.reserve(coordinates.size());
        for (size_t n = 0; n < coordinates.size(); n++)
        {
            target.push_back(n);
        }
    }
    else
    {
        std::vector<size_t> indices;
        indices.reserve(coordinates.size());
        for (auto& a : coordinates)
        {
            indices.push_back(gridIndex(bounds, a));
        }
        auto& intermediate = this->intermediate_;
        intermediate.resize(size);
        for (auto& a : intermediate)
        {
            a.clear();
        }
        for (size_t n = 0; n < indices.size(); n++)
        {
            intermediate[indices[n]].push_back(n);
        }

        auto insert = [](std::vector<size_t>& target, std::vector<size_t>& origin)
        {
            target.insert(target.end(), origin.begin(), origin.end());
            return true;
        };

        for (size_t n = 0; n < size; n++)
        {
            result[n].clear();
            auto& target = result[n];
            for (auto offset : gridNeighborOffsets(subdivisions, bounds.stride))
            {
                size_t k =
                    n + offset; // can underflow at negative offsets, and that's fine since (k < size) will be false
                (k < size) && insert(target, intermediate[k]);
            }
        }
    }
}
