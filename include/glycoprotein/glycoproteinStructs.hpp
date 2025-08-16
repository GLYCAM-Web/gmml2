#ifndef INCLUDE_GLYCOPROTEIN_GLYCOPROTEINSTRUCTS_HPP
#define INCLUDE_GLYCOPROTEIN_GLYCOPROTEINSTRUCTS_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralAngleSearchTypes.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/geometry/geometryTypes.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/metadata/residueTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        enum class MoleculeType
        {
            protein,
            glycan
        };

        struct AtomData
        {
            std::vector<std::string> names;
            std::vector<std::string> types;
            std::vector<uint> numbers;
            std::vector<uint> serializedNumbers;
            std::vector<uint> atomicNumbers;
            std::vector<std::string> elementStrings;
            std::vector<Element> elements;
            std::vector<double> charges;
            std::vector<Sphere> initialState;
            std::vector<bool> all;
            std::vector<bool> includeInEachOverlapCheck;
            std::vector<bool> includeInMainOverlapCheck;
            std::vector<bool> partOfMovableSidechain;
        };

        struct SidechainDihedral
        {
            std::array<size_t, 4> atoms;
            std::vector<size_t> movingAtoms;
        };

        struct ResidueData
        {
            std::vector<std::string> names;
            std::vector<ResidueType> types;
            std::vector<bool> hasAllExpectedAtoms;
            std::vector<double> phiAngles;
            std::vector<double> psiAngles;
            std::vector<std::string> ids;
            std::vector<uint> numbers;
            std::vector<uint> serializedNumbers;
            std::vector<std::string> chainIds;
            std::vector<std::vector<SidechainDihedral>> sidechainDihedrals;
            std::vector<std::vector<size_t>> sidechainRotations;
            std::vector<std::vector<double>> sidechainRotationWeights;
            std::vector<Sphere> sidechainPotentialBounds;
        };

        struct MoleculeData
        {
            std::vector<MoleculeType> types;
            std::vector<bool> isProtein;
        };

        struct ResidueEdgeData
        {
            std::vector<std::array<std::vector<bool>, 2>> atomsCloseToEdge;
        };

        struct RotatableBondData
        {
            std::vector<std::array<size_t, 4>> dihedralAtoms;
            std::vector<std::vector<size_t>> movingAtoms;
            std::vector<std::vector<size_t>> metadata;
            std::vector<AngleWithMetadata> initialShape;
        };

        struct ResidueLinkageData
        {
            std::vector<RotamerType> rotamerTypes;
            std::vector<bool> isGlycositeLinkage;
            std::vector<size_t> edgeId;
            std::vector<std::vector<size_t>> rotatableBonds;
        };

        struct GlycanData
        {
            std::vector<size_t> attachmentResidue;
            std::vector<size_t> moleculeId;
            std::vector<std::vector<size_t>> linkages;
        };

        struct AssemblyData
        {
            AtomData atoms;
            ResidueData residues;
            MoleculeData molecules;
            ResidueEdgeData residueEdges;
            RotatableBondData rotatableBonds;
            ResidueLinkageData residueLinkages;
            GlycanData glycans;
        };

        struct MutableData
        {
            assembly::Bounds bounds;
            std::vector<size_t> rotatableBondCurrentMetadata;
            std::vector<bool> moleculeIncluded;
            std::vector<bool> residueSidechainMoved;
        };

        struct OverlapSettings
        {
            PotentialTable potentialTable;
            double tolerance;
            double rejectionThreshold;
        };

        struct GlycoproteinAssembly
        {
            assembly::Graph graph;
            AssemblyData data;
            MutableData mutableData;
        };
    } // namespace gpbuilder
} // namespace gmml

#endif
