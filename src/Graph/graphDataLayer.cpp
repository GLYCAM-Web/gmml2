#include "includes/Graph/graphDataLayer.hpp"
#include "includes/Graph/manipulation.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"

#include <cstddef>
#include <optional>
#include <vector>
#include <string>
#include <stdexcept>

namespace graph
{
    size_t GraphDataLayer::addAtom(size_t residue, cds::Coordinate coord, std::string name, std::string element,
                                   double charge)
    {
        atomResidue.push_back(residue);
        atomCoordinates.push_back(coord);
        atomNames.push_back(name);
        atomElements.push_back(element);
        atomCharges.push_back(charge);
        return addNode(database);
    }

    size_t GraphDataLayer::addResidue(size_t molecule, std::string name)
    {
        size_t index = residueMolecule.size();
        residueMolecule.push_back(molecule);
        residueNames.push_back(name);
        return index;
    }

    size_t GraphDataLayer::addMolecule(std::string name)
    {
        size_t index = moleculeNames.size();
        moleculeNames.push_back(name);
        return index;
    }

    size_t GraphDataLayer::addBond(size_t sourceAtom, size_t targetAtom, BondType type,
                                   std::optional<ResidueLinkageStruct> residueLinkage,
                                   std::optional<MoleculeLinkageStruct> moleculeLinkage)
    {
        size_t sourceResidue    = atomResidue[sourceAtom];
        size_t targetResidue    = atomResidue[targetAtom];
        bool differentResidues  = sourceResidue != targetResidue;
        bool differentMolecules = residueMolecule[sourceResidue] != residueMolecule[targetResidue];
        if (residueLinkage.has_value() != differentResidues)
        {
            if (differentResidues)
            {
                throw std::runtime_error("Missing residue linkage for bond between different residues");
            }
            else
            {
                throw std::runtime_error("Given residue linkage for bond within same residue");
            }
        }
        if (moleculeLinkage.has_value() != differentMolecules)
        {
            if (differentMolecules)
            {
                throw std::runtime_error("Missing molecule linkage for bond between different molecules");
            }
            else
            {
                throw std::runtime_error("Given molecule linkage for bond within same molecule");
            }
        }
        bondTypes.push_back(type);
        residueLinkages.push_back(residueLinkage);
        moleculeLinkages.push_back(moleculeLinkage);
        return addEdge(database, {sourceAtom, targetAtom});
    }

    void GraphDataLayer::removeBond(size_t index)
    {
        removeEdge(database, index);
    }

    void GraphDataLayer::removeAtom(size_t index)
    {
        removeNode(database, index);
    }

    void GraphDataLayer::removeResidue(size_t index)
    {
        for (size_t n = 0; n < atomResidue.size(); n++)
        {
            if (atomResidue[n] == index)
            {
                removeAtom(n);
            }
        }
    }

    void GraphDataLayer::removeMolecule(size_t index)
    {
        for (size_t n = 0; n < residueMolecule.size(); n++)
        {
            if (residueMolecule[n] == index)
            {
                removeResidue(n);
            }
        }
    }
} // namespace graph
