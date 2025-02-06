#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SIDECHAINS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SIDECHAINS_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/MolecularMetadata/sidechainRotamers.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    bool sidechainHasGlycanOverlap(const assembly::Graph& graph, const AssemblyData& data,
                                   const MutableData& mutableData, size_t sidechainResidue);
    void updateSidechainRotation(const MolecularMetadata::SidechainRotamerData& sidechains,
                                 const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                                 size_t sidechainResidue, size_t rotation);
    size_t lowestOverlapSidechainRotation(const MolecularMetadata::SidechainRotamerData& sidechains,
                                          const assembly::Graph& graph, const AssemblyData& data,
                                          MutableData& mutableData, size_t sidechainResidue,
                                          const std::vector<size_t>& otherAtoms);
    std::vector<size_t> atomsWithinSidechainPotentialBounds(const assembly::Graph& graph, const AssemblyData& data,
                                                            const MutableData& mutableData, size_t sidechainResidue);
    std::vector<std::vector<SidechainDihedral>> sidechainDihedrals(const assembly::Graph& graph,
                                                                   const AssemblyData& data);
    std::vector<std::vector<size_t>> sidechainRotations(const assembly::Graph& graph, const AssemblyData& data,
                                                        const MolecularMetadata::SidechainRotamerData& sidechains);
    std::vector<cds::Sphere> sidechainPotentialBounds(const assembly::Graph& graph, const AssemblyData& data,
                                                      const MutableData& mutableData,
                                                      const MolecularMetadata::SidechainRotamerData& sidechains);
    GlycoproteinAssembly addSidechainRotamers(const MolecularMetadata::SidechainRotamerData& sidechains,
                                              GlycoproteinAssembly assembly);
    std::vector<bool> partOfMovableSidechain(const assembly::Graph& graph, const AssemblyData& data);
} // namespace glycoproteinBuilder
#endif
