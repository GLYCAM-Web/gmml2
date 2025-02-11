#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP

#include "includes/Graph/graphTypes.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    enum class MoleculeType
    {
        protein,
        glycan
    };

    struct OverlapWeight
    {
        double protein;
        double glycan;
        double self;
    };

    struct AtomData
    {
        std::vector<std::string> names;
        std::vector<std::string> types;
        std::vector<int> numbers;
        std::vector<int> serializedNumbers;
        std::vector<int> atomicNumbers;
        std::vector<std::string> elements;
        std::vector<double> charges;
        std::vector<cds::Sphere> initialState;
        std::vector<bool> partOfMovableSidechain;
        std::vector<bool> all;
        std::vector<bool> alwaysIncluded;
    };

    struct SidechainDihedral
    {
        std::array<size_t, 4> atoms;
        std::vector<size_t> movingAtoms;
    };

    struct ResidueData
    {
        std::vector<std::string> names;
        std::vector<cds::ResidueType> types;
        std::vector<bool> hasAllExpectedAtoms;
        std::vector<double> phiAngles;
        std::vector<double> psiAngles;
        std::vector<std::string> ids;
        std::vector<int> numbers;
        std::vector<int> serializedNumbers;
        std::vector<double> overlapWeights;
        std::vector<std::vector<SidechainDihedral>> sidechainDihedrals;
        std::vector<std::vector<size_t>> sidechainRotations;
        std::vector<cds::Sphere> sidechainPotentialBounds;
    };

    struct MoleculeData
    {
        std::vector<MoleculeType> types;
    };

    struct GlycanData
    {
        std::vector<size_t> attachmentResidue;
        std::vector<size_t> moleculeId;
        std::vector<std::vector<size_t>> linkages;
    };

    struct RotatableDihedralData
    {
        std::vector<GlycamMetadata::DihedralAngleDataVector> metadata;
        std::vector<cds::AngleWithMetadata> initialShape;
    };

    struct ResidueLinkageData
    {
        std::vector<GlycamMetadata::RotamerType> rotamerTypes;
        std::vector<std::vector<cds::BondedResidueOverlapInput>> overlapBonds;
        std::vector<bool> branching;
        std::vector<bool> isGlycositeLinkage;
    };

    struct RotatableDihedralIndices
    {
        std::array<size_t, 4> atoms;
        std::vector<size_t> movingAtoms;
    };

    struct ResidueLinkageIndices
    {
        size_t residueEdge;
        std::vector<size_t> rotatableDihedrals;
        std::array<std::vector<bool>, 2> closelyBondedAtoms;
        std::vector<size_t> nonReducingResidues;
        std::vector<size_t> reducingResidues;
    };

    struct AssemblyIndices
    {
        std::vector<size_t> proteinMolecules;
        std::vector<RotatableDihedralIndices> rotatableDihedrals;
        std::vector<ResidueLinkageIndices> residueLinkages;
    };

    struct AssemblyData
    {
        AtomData atoms;
        ResidueData residues;
        MoleculeData molecules;
        GlycanData glycans;
        RotatableDihedralData rotatableDihedralData;
        ResidueLinkageData residueLinkageData;
        AssemblyIndices indices;
    };

    struct MutableData
    {
        std::vector<cds::Sphere> atomBounds;
        std::vector<cds::Sphere> residueBounds;
        std::vector<cds::Sphere> moleculeBounds;
        std::vector<bool> glycanIncluded;
        std::vector<bool> residueSidechainMoved;
    };

    struct GlycoproteinAssembly
    {
        assembly::Graph graph;
        AssemblyData data;
        MutableData mutableData;
    };

    typedef std::function<void(pcg32&, const assembly::Graph&, const AssemblyData&, MutableData&,
                               const std::vector<std::vector<cds::ResidueLinkageShapePreference>>&,
                               const std::vector<size_t>&)>
        SidechainAdjustment;

} // namespace glycoproteinBuilder
#endif
