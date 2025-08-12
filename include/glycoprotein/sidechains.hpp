#ifndef INCLUDE_GLYCOPROTEIN_SIDECHAINS_HPP
#define INCLUDE_GLYCOPROTEIN_SIDECHAINS_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/metadata/aminoAcids.hpp"
#include "include/metadata/sidechainRotamers.hpp"

#include <vector>

namespace gmml
{
    namespace gpbuilder
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
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const MutableData& mutableData,
            const std::vector<size_t>& glycans,
            size_t sidechainResidue);

        void updateSidechainRotation(
            const SidechainRotamerData& sidechains,
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t sidechainResidue,
            size_t rotation);

        void restoreSidechainRotation(
            const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData, size_t residue);

        void setSidechainToLowestOverlapState(
            const SidechainRotamerData& sidechains,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            const std::vector<size_t>& preference,
            size_t residue);

        std::vector<double> sidechainOverlap(
            const OverlapSettings& overlapSettings,
            const assembly::Indices& indices,
            const AssemblyData& data,
            const std::vector<Sphere>& bounds,
            const std::vector<size_t>& atomsA,
            const std::vector<size_t>& atomsB);

        IndexedOverlap lowestOverlapSidechainRotation(
            const SidechainRotamerData& sidechains,
            const OverlapSettings& overlapSettings,
            const assembly::Indices& indices,
            const AssemblyData& data,
            const MutableData& mutableData,
            const std::vector<size_t>& preference,
            size_t sidechainResidue,
            const std::vector<size_t>& otherAtoms);

        std::vector<size_t> atomsWithinSidechainPotentialBounds(
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const MutableData& mutableData,
            const std::vector<bool>& includedAtoms,
            size_t sidechainResidue);

        std::vector<std::vector<SidechainDihedral>> sidechainDihedrals(
            const AminoAcidTable& aminoAcidTable, const assembly::Graph& graph, const AssemblyData& data);

        SidechainRotationsAndWeights sidechainRotationsAndWeights(
            const assembly::Indices& indices, const AssemblyData& data, const SidechainRotamerData& sidechains);

        std::vector<Sphere> sidechainPotentialBounds(
            const assembly::Indices& indices,
            const AssemblyData& data,
            const MutableData& mutableData,
            const SidechainRotamerData& sidechains);

        GlycoproteinAssembly addSidechainRotamers(
            const AminoAcidTable& aminoAcidTable,
            const SidechainRotamerData& sidechains,
            GlycoproteinAssembly assembly);

        std::vector<bool> partOfMovableSidechain(const assembly::Indices& indices, const AssemblyData& data);
    } // namespace gpbuilder
} // namespace gmml

#endif
