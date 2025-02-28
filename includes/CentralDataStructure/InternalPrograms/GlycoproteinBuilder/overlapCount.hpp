#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblySelection.hpp"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    std::vector<cds::Overlap> intraGlycanOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                                  const assembly::Selection& selection, const assembly::Bounds& bounds,
                                                  const std::vector<double>& residueWeights, size_t glycanId);
    std::vector<cds::Overlap> moleculeOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                               const assembly::Selection& selection, const assembly::Bounds& bounds,
                                               const std::vector<double>& residueWeights, size_t moleculeA,
                                               size_t moleculeB);
    std::vector<cds::Overlap> moleculeResidueOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                                      const assembly::Selection& selection,
                                                      const assembly::Bounds& bounds,
                                                      const std::vector<double>& residueWeights, size_t molecule,
                                                      size_t residue);
    std::vector<cds::Overlap> totalOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                            const assembly::Selection& selection, const assembly::Bounds& bounds,
                                            const std::vector<double>& residueWeights, OverlapMultiplier weight);
    cds::Overlap localOverlap(const assembly::Graph& graph, const AssemblyData& data,
                              const assembly::Selection& selection, const assembly::Bounds& bounds,
                              const std::vector<double>& residueWeights, size_t glycanId, double selfWeight);
    std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const assembly::Graph& graph,
                                                  const AssemblyData& data, const assembly::Selection& selection,
                                                  const assembly::Bounds& bounds);
} // namespace glycoproteinBuilder
#endif
