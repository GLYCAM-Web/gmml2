#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SIDECHAINS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SIDECHAINS_HPP

#include "includes/Assembly/assemblyTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/MolecularMetadata/sidechainRotamers.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    struct IndexedOverlap
    {
        size_t index;
        double overlap;
    };

    struct SidechainRotationsAndWeights
    {
        std::vector<std::vector<size_t>> rotations;
        std::vector<std::vector<double>> weights;
    };

    bool sidechainHasGlycanOverlap(
        const assembly::Graph& graph,
        const AssemblyData& data,
        const MutableData& mutableData,
        const std::vector<size_t>& glycans,
        size_t sidechainResidue);

    void updateSidechainRotation(
        const MolecularMetadata::SidechainRotamerData& sidechains,
        const assembly::Graph& graph,
        const AssemblyData& data,
        MutableData& mutableData,
        size_t sidechainResidue,
        size_t rotation);

    void restoreSidechainRotation(
        const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData, size_t residue);

    void setSidechainToLowestOverlapState(
        const MolecularMetadata::SidechainRotamerData& sidechains,
        const assembly::Graph& graph,
        const AssemblyData& data,
        MutableData& mutableData,
        const std::vector<size_t>& preference,
        size_t residue);

    std::vector<double> sidechainOverlap(
        const assembly::Indices& indices,
        const AssemblyData& data,
        const std::vector<cds::Sphere>& bounds,
        const std::vector<size_t>& atomsA,
        const std::vector<size_t>& atomsB);

    IndexedOverlap lowestOverlapSidechainRotation(
        const MolecularMetadata::SidechainRotamerData& sidechains,
        const assembly::Indices& indices,
        const AssemblyData& data,
        const MutableData& mutableData,
        const std::vector<size_t>& preference,
        size_t sidechainResidue,
        const std::vector<size_t>& otherAtoms);

    std::vector<size_t> atomsWithinSidechainPotentialBounds(
        const assembly::Graph& graph,
        const AssemblyData& data,
        const MutableData& mutableData,
        const std::vector<bool>& includedAtoms,
        size_t sidechainResidue);

    std::vector<std::vector<SidechainDihedral>> sidechainDihedrals(
        const MolecularMetadata::AminoAcidTable& aminoAcidTable,
        const assembly::Graph& graph,
        const AssemblyData& data);

    SidechainRotationsAndWeights sidechainRotationsAndWeights(
        const assembly::Indices& indices,
        const AssemblyData& data,
        const MolecularMetadata::SidechainRotamerData& sidechains);

    std::vector<cds::Sphere> sidechainPotentialBounds(
        const assembly::Indices& indices,
        const AssemblyData& data,
        const MutableData& mutableData,
        const MolecularMetadata::SidechainRotamerData& sidechains);

    GlycoproteinAssembly addSidechainRotamers(
        const MolecularMetadata::AminoAcidTable& aminoAcidTable,
        const MolecularMetadata::SidechainRotamerData& sidechains,
        GlycoproteinAssembly assembly);

    std::vector<bool> partOfMovableSidechain(const assembly::Indices& indices, const AssemblyData& data);
} // namespace glycoproteinBuilder
#endif
