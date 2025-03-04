#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP

#include "includes/Graph/graphTypes.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

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
        std::vector<bool> partOfMovableSidechain;
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
    };

    struct MoleculeData
    {
        std::vector<MoleculeType> types;
    };

    struct RotatableDihedralData
    {
        std::vector<GlycamMetadata::DihedralAngleDataVector> metadata;
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

    struct GlycanIndices
    {
        size_t attachmentResidue;
        size_t glycanMolecule;
        std::vector<size_t> linkages;
    };

    struct AssemblyIndices
    {
        std::vector<size_t> proteinMolecules;
        std::vector<RotatableDihedralIndices> rotatableDihedrals;
        std::vector<ResidueLinkageIndices> residueLinkages;
        std::vector<GlycanIndices> glycans;
    };

    struct AssemblyData
    {
        AtomData atoms;
        ResidueData residues;
        MoleculeData molecules;
        RotatableDihedralData rotatableDihedralData;
        ResidueLinkageData residueLinkageData;
        AssemblyIndices indices;
    };

    struct MutableData
    {
        std::vector<cds::Sphere> atomBounds;
        std::vector<cds::Sphere> residueBounds;
        std::vector<cds::Sphere> moleculeBounds;
        std::vector<cds::AngleWithMetadata> currentDihedralShape;
        std::vector<bool> atomIgnored;
        std::vector<bool> glycanIncluded;
    };

    struct GlycoproteinAssembly
    {
        assembly::Graph graph;
        AssemblyData data;
        MutableData mutableData;
    };

} // namespace glycoproteinBuilder
#endif
