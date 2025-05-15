#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBDATA_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBDATA_HPP

#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/Graph/graphTypes.hpp"

#include <string>
#include <vector>

namespace pdb
{
    struct AtomEntry
    {
        std::string recordName;
        std::string name;
        uint number;
        cds::Coordinate coordinate;
        double occupancy;
        double temperatureFactor;
    };

    struct AtomData
    {
        std::vector<std::string> recordNames;
        std::vector<uint> numbers;
        std::vector<std::string> names;
        std::vector<MolecularMetadata::Element> elements;
        std::vector<cds::Coordinate> coordinates;
        std::vector<double> occupancies;
        std::vector<double> temperatureFactors;
    };

    struct ResidueEntry
    {
        std::string name;
        cds::ResidueType type;
        uint number;
        std::string insertionCode;
        std::string chainId;
        bool hasTerCard;
    };

    struct ResidueData
    {
        std::vector<std::string> names;
        std::vector<cds::ResidueType> types;
        std::vector<uint> numbers;
        std::vector<std::string> insertionCodes;
        std::vector<std::string> chainIds;
        std::vector<bool> isCTerminal;
        std::vector<bool> isNTerminal;
        std::vector<bool> hasTerCard;
    };

    struct MoleculeData
    {
        std::vector<std::string> chainIds;
        std::vector<std::vector<size_t>> residueOrder;
    };

    struct AssemblyData
    {
        std::vector<uint> numbers;
    };

    struct PdbData
    {
        AtomData atoms;
        ResidueData residues;
        MoleculeData molecules;
        AssemblyData assemblies;
        graph::Database atomGraph;
        assembly::Indices indices;
        cds::GraphObjects objects;
    };
} // namespace pdb

#endif