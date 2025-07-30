#ifndef INCLUDE_PDB_PDBDATA_HPP
#define INCLUDE_PDB_PDBDATA_HPP

#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/metadata/residueTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        enum ResidueTerminality
        {
            NonTerminal = 0,
            NTerminal = 1,
            CTerminal = 2
        };

        static const std::vector<std::string> terminalPrefix {"", "N", "C"};

        struct AtomEntry
        {
            std::string recordName;
            std::string name;
            uint number;
            Coordinate coordinate;
            double occupancy;
            double temperatureFactor;
        };

        struct AtomData
        {
            std::vector<std::string> recordNames;
            std::vector<uint> numbers;
            std::vector<std::string> names;
            std::vector<Element> elements;
            std::vector<Coordinate> coordinates;
            std::vector<double> occupancies;
            std::vector<double> temperatureFactors;
            std::vector<double> charges;
            std::vector<std::string> types;
        };

        struct ResidueEntry
        {
            std::string name;
            ResidueType type;
            uint number;
            std::string insertionCode;
            std::string chainId;
            bool hasTerCard;
        };

        struct ResidueData
        {
            std::vector<std::string> names;
            std::vector<ResidueType> types;
            std::vector<uint> numbers;
            std::vector<std::string> insertionCodes;
            std::vector<std::string> chainIds;
            std::vector<ResidueTerminality> terminality;
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

        struct TrajectoryData
        {
            std::vector<std::vector<Coordinate>> coordinates;
        };

        struct PdbData
        {
            AtomData atoms;
            ResidueData residues;
            MoleculeData molecules;
            AssemblyData assemblies;
            TrajectoryData trajectory;
            graph::Database atomGraph;
            assembly::Indices indices;
            GraphObjects objects;
        };
    } // namespace pdb
} // namespace gmml

#endif