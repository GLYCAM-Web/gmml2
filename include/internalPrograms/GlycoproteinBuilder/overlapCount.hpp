#ifndef INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP
#define INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"

#include <functional>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        std::vector<double> totalOverlaps(
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds);

        double localOverlap(
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds,
            size_t glycanId);

        std::vector<size_t> determineSitesWithOverlap(
            double threshold,
            const std::vector<size_t>& movedSites,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds);
    } // namespace gpbuilder
} // namespace gmml

#endif
