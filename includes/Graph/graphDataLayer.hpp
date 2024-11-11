#ifndef INCLUDES_CENTRALDATASTRUCTURE_GRAPH_GRAPHDATALAYER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GRAPH_GRAPHDATALAYER_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/Graph/types.hpp"

#include <cstddef>
#include <array>
#include <optional>
#include <vector>
#include <string>

namespace graph
{
    enum class BondType
    {
        covalent,
        james
    };

    struct ResidueLinkageStruct
    {
        std::vector<size_t> rotatableDihedrals;
    };

    struct MoleculeLinkageStruct
    {};

    class GraphDataLayer
    {
      public:
        size_t addAtom(size_t residue, cds::Coordinate coord, std::string name, std::string element, double charge);
        size_t addResidue(size_t molecule, std::string name);
        size_t addMolecule(std::string name);
        size_t addBond(size_t sourceAtom, size_t targetAtom, BondType type,
                       std::optional<ResidueLinkageStruct> residueLinkage,
                       std::optional<MoleculeLinkageStruct> moleculeLinkage);
        void removeBond(size_t index);
        void removeAtom(size_t index);
        void removeResidue(size_t index);
        void removeMolecule(size_t index);

        graph::Database database;

        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        std::vector<cds::Coordinate> atomCoordinates;
        std::vector<std::string> atomNames;
        std::vector<std::string> atomElements;
        std::vector<double> atomCharges;
        std::vector<std::string> residueNames;
        std::vector<std::string> moleculeNames;
        std::vector<BondType> bondTypes;
        std::vector<std::optional<ResidueLinkageStruct>> residueLinkages;
        std::vector<std::optional<MoleculeLinkageStruct>> moleculeLinkages;
    };
} // namespace graph

#endif
