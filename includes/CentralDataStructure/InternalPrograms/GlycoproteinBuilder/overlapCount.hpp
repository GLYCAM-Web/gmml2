#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    std::vector<cds::Overlap> intraGlycanOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                                  const MutableData& mutableData,
                                                  const std::vector<double>& residueWeights,
                                                  const std::vector<bool>& includedAtoms, size_t glycanId);
    std::vector<cds::Overlap> moleculeOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                               const MutableData& mutableData,
                                               const std::vector<double>& residueWeights,
                                               const std::vector<bool>& includedAtoms, size_t moleculeA,
                                               size_t moleculeB);
    std::vector<cds::Overlap> moleculeResidueOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                                      const MutableData& mutableData,
                                                      const std::vector<double>& residueWeights,
                                                      const std::vector<bool>& includedAtoms, size_t molecule,
                                                      size_t residue);
    std::vector<cds::Overlap> totalOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                            const MutableData& mutableData, const std::vector<double>& residueWeights,
                                            const std::vector<bool>& includedAtoms, OverlapMultiplier weight);
    cds::Overlap localOverlap(const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData,
                              const std::vector<double>& residueWeights, const std::vector<bool>& includedAtoms,
                              size_t glycanId, double selfWeight);
    std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const assembly::Graph& graph,
                                                  const AssemblyData& data, const MutableData& mutableData,
                                                  const std::vector<bool>& includedAtoms);
} // namespace glycoproteinBuilder
#endif
