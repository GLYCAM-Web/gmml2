#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_BINARYPARTITIONGRID_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_BINARYPARTITIONGRID_HPP

#include "includes/CentralDataStructure/coordinate.hpp"

#include <array>
#include <vector>

namespace cds
{
    struct BinaryPartitionGridAxisLimits
    {
        // acts as a multidimensional array
        // first value is from first subdivision, next two from second, etc
        std::vector<double> values;
    };

    typedef std::array<BinaryPartitionGridAxisLimits, 3> GridBoundsAxes;
    typedef std::array<size_t, 3> GridStride;

    struct BinaryPartitionGridBounds
    {
        size_t subdivisions;
        GridStride stride;
        GridBoundsAxes limits;
    };

    struct BinaryPartitionGrid
    {
        BinaryPartitionGridBounds bounds;
        std::vector<std::vector<size_t>> intermediate;
        std::vector<std::vector<size_t>> grid;
    };

    inline size_t gridIndexPosition(size_t subdivisions, const BinaryPartitionGridAxisLimits& limits, double x)
    {
        size_t index  = 0;
        size_t offset = 0;
        size_t width  = 1;
        for (size_t n = 0; n < subdivisions; n++)
        {
            index  = (2 * index) + (x > limits.values[index + offset]);
            offset += width;
            width  *= 2;
        }
        return index;
    }

    inline size_t gridIndex(const BinaryPartitionGridBounds& boundary, const cds::Coordinate& coord)
    {
        auto at = [&](size_t n)
        {
            size_t position = gridIndexPosition(boundary.subdivisions, boundary.limits[n], coord.nth(n));
            return boundary.stride[n] * position;
        };
        return at(0) + at(1) + at(2);
    }

    inline const std::vector<size_t>& gridAt(const BinaryPartitionGrid& grid, const cds::Coordinate& coord)
    {
        return grid.grid[gridIndex(grid.bounds, coord)];
    }

    BinaryPartitionGridBounds binaryPartitionGridBounds(size_t subdivisions, double distance, const Coordinate& coord);
    void initBinaryPartitionGrid(BinaryPartitionGrid& grid, const BinaryPartitionGridBounds bounds,
                                 const std::vector<Coordinate>& coordinates);
} // namespace cds
#endif
