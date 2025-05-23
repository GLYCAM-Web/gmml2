#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/Assembly/assemblyTypes.hpp"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    std::vector<double> totalOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                      const assembly::Selection& selection, const assembly::Bounds& bounds);
    double localOverlap(const assembly::Graph& graph, const AssemblyData& data, const assembly::Selection& selection,
                        const assembly::Bounds& bounds, size_t glycanId);
    std::vector<size_t> determineSitesWithOverlap(double threshold, const std::vector<size_t>& movedSites,
                                                  const assembly::Graph& graph, const AssemblyData& data,
                                                  const assembly::Selection& selection, const assembly::Bounds& bounds);
} // namespace glycoproteinBuilder
#endif
