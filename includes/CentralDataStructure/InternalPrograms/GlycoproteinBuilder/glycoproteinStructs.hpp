#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP

#include "includes/Graph/graphTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
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
        std::vector<int> atomicNumbers;
        std::vector<std::string> elements;
        std::vector<double> charges;
        std::vector<cds::Sphere> bounds;
    };

    struct ResidueData
    {
        std::vector<std::string> names;
        std::vector<cds::ResidueType> types;
        std::vector<int> numbers;
        std::vector<double> overlapWeights;
        std::vector<cds::Sphere> bounds;
    };

    struct MoleculeData
    {
        std::vector<MoleculeType> types;
        std::vector<cds::Sphere> bounds;
    };

    struct RotatableDihedralData
    {
        std::vector<cds::AngleWithMetadata> currentShape;
    };

    struct ResidueLinkagedata
    {
        std::vector<GlycamMetadata::RotamerType> rotamerTypes;
        std::vector<std::vector<GlycamMetadata::DihedralAngleDataVector>> metadata;
        std::vector<std::vector<cds::BondedResidueOverlapInput>> overlapBonds;
        std::vector<bool> branching;
    };

    struct GlycanData
    {
        std::vector<bool> included;
    };

    struct AssemblyData
    {
        AtomData atoms;
        ResidueData residues;
        MoleculeData molecules;
        RotatableDihedralData rotatableDihedralData;
        ResidueLinkagedata residueLinkageData;
        GlycanData glycanData;
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

    struct AssemblyGraphs
    {
        cds::GraphIndexData indices;
        graph::Graph atoms;
        graph::Graph residues;
        graph::Graph molecules;
        std::vector<size_t> proteinMolecules;
        std::vector<RotatableDihedralIndices> rotatableDihedralIndices;
        std::vector<ResidueLinkageIndices> residueLinkages;
        std::vector<GlycanIndices> glycans;
    };

    struct GlycoproteinAssembly
    {
        AssemblyGraphs graphs;
        AssemblyData data;
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
