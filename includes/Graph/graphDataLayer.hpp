#ifndef INCLUDES_CENTRALDATASTRUCTURE_GRAPH_GRAPHDATALAYER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GRAPH_GRAPHDATALAYER_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/Graph/types.hpp"
#include "includes/Graph/manipulation.hpp"

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
    {};

    struct MoleculeLinkageStruct
    {};

    class GraphDataLayer
    {
      public:
        inline size_t addAtom(size_t residue, cds::Coordinate coord, std::string name, std::string element,
                              double charge)
        {
            atomResidue.push_back(residue);
            atomCoordinates.push_back(coord);
            atomNames.push_back(name);
            atomElements.push_back(element);
            atomCharges.push_back(charge);
            return addNode(database);
        }

        inline size_t addResidue(size_t molecule, std::string name)
        {
            size_t index = residueMolecule.size();
            residueMolecule.push_back(molecule);
            residueNames.push_back(name);
            return index;
        }

        inline size_t addMolecule(std::string name)
        {
            size_t index = moleculeNames.size();
            moleculeNames.push_back(name);
            return index;
        }

        inline size_t addBond(size_t sourceAtom, size_t targetAtom, BondType type,
                              std::optional<ResidueLinkageStruct> residueLinkage,
                              std::optional<MoleculeLinkageStruct> moleculeLinkage)
        {
            // assertions
            // size_t sourceResidue = atomResidue[sourceAtom];
            // size_t targetResidue = atomResidue[targetAtom];
            // residueLinkage.has_value() == (sourceResidue != targetResidue)
            // moleculeLinkage.has_value() == (residueMolecule[sourceResidue] != residueMolecule[targetResidue])
            bondTypes.push_back(type);
            residueLinkages.push_back(residueLinkage);
            moleculeLinkages.push_back(moleculeLinkage);
            return addEdge(database, {sourceAtom, targetAtom});
        }

        inline void removeBond(size_t index)
        {
            removeEdge(database, index);
        }

        inline void removeAtom(size_t index)
        {
            removeNode(database, index);
        }

        inline void removeResidue(size_t index)
        {
            for (size_t n = 0; n < atomResidue.size(); n++)
            {
                if (atomResidue[n] == index)
                {
                    removeAtom(n);
                }
            }
        }

        inline void removeMolecule(size_t index)
        {
            for (size_t n = 0; n < residueMolecule.size(); n++)
            {
                if (residueMolecule[n] == index)
                {
                    removeResidue(n);
                }
            }
        }

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
