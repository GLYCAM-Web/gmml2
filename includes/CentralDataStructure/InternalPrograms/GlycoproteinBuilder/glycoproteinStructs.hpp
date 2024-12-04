#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP

#include "includes/Graph/graphTypes.hpp"
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
    };

    struct ResidueData
    {
        std::vector<std::string> names;
        std::vector<cds::ResidueType> types;
        std::vector<std::string> ids;
        std::vector<int> numbers;
        std::vector<int> serializedNumbers;
        std::vector<double> overlapWeights;
    };

    struct MoleculeData
    {
        std::vector<MoleculeType> types;
    };

    struct ResidueLinkageData
    {
        std::vector<GlycamMetadata::RotamerType> rotamerTypes;
        std::vector<std::vector<GlycamMetadata::DihedralAngleDataVector>> metadata;
        std::vector<std::vector<cds::BondedResidueOverlapInput>> overlapBonds;
        std::vector<bool> branching;
        std::vector<bool> isGlycositeLinkage;
    };

    struct AssemblyData
    {
        AtomData atoms;
        ResidueData residues;
        MoleculeData molecules;
        ResidueLinkageData residueLinkageData;
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
        size_t atomCount;
        size_t residueCount;
        size_t moleculeCount;
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        std::vector<size_t> proteinMolecules;
        std::vector<RotatableDihedralIndices> rotatableDihedrals;
        std::vector<ResidueLinkageIndices> residueLinkages;
        std::vector<GlycanIndices> glycans;
    };

    struct AssemblyGraphs
    {
        AssemblyIndices indices;
        graph::Graph atoms;
        graph::Graph residues;
        graph::Graph molecules;
    };

    struct MutableData
    {
        std::vector<cds::Sphere> atomBounds;
        std::vector<cds::Sphere> residueBounds;
        std::vector<cds::Sphere> moleculeBounds;
        std::vector<cds::AngleWithMetadata> currentDihedralShape;
        std::vector<bool> glycanIncluded;
    };

    struct GlycoproteinAssembly
    {
        AssemblyGraphs graphs;
        AssemblyData data;
        MutableData mutableData;
    };

    inline const std::vector<size_t>& residueAtoms(const AssemblyGraphs& graphs, size_t residueId)
    {
        return graphs.residues.nodes.elements[residueId];
    };

    inline const std::vector<size_t>& moleculeResidues(const AssemblyGraphs& graphs, size_t moleculeId)
    {
        return graphs.molecules.nodes.elements[moleculeId];
    };

    inline size_t moleculeEdgeToAtomEdgeIndex(const AssemblyGraphs& graphs, size_t moleculeEdgeId)
    {
        size_t residueBondIndex = graphs.molecules.edges.indices[moleculeEdgeId];
        return graphs.residues.edges.indices[residueBondIndex];
    }

} // namespace glycoproteinBuilder
#endif
