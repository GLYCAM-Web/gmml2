#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPCOUNT_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    cds::Overlap intraGlycanOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                     const MutableData& mutableData, size_t glycanId);
    cds::Overlap moleculeOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                  const MutableData& mutableData, size_t moleculeA, size_t moleculeB);
    cds::Overlap moleculeResidueOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                         const MutableData& mutableData, size_t molecule, size_t residue);
    cds::Overlap totalOverlaps(const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData,
                               OverlapWeight weight);
    cds::Overlap localOverlap(const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData,
                              size_t glycanId, double selfWeight);
    std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const assembly::Graph& graph,
                                                  const AssemblyData& data, const MutableData& mutableData);
} // namespace glycoproteinBuilder
#endif
