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

        std::vector<size_t> determineSitesWithOverlap(
            double threshold,
            const OverlapSettings& overlapSettings,
            const std::vector<size_t>& movedSites,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds);
    } // namespace gpbuilder
} // namespace gmml

#endif
