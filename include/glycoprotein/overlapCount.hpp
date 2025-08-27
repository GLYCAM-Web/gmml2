#ifndef INCLUDE_GLYCOPROTEIN_OVERLAPCOUNT_HPP
#define INCLUDE_GLYCOPROTEIN_OVERLAPCOUNT_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"

#include <functional>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        struct OverlapSites
        {
            std::vector<size_t> indices;
            std::vector<bool> aboveOverlapThreshold;
            std::vector<std::vector<bool>> interactions;
            std::vector<size_t> concertIds;
        };

        std::vector<double> totalOverlaps(
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds);

        double localOverlap(
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds,
            size_t glycanId);

        OverlapSites determineOverlapState(
            double threshold,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds);
    } // namespace gpbuilder
} // namespace gmml

#endif
