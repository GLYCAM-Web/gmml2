#include "includes/CentralDataStructure/Selections/residueSelections.hpp"

#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/containers.hpp"

using cds::Residue;

std::vector<Residue*> cdsSelections::selectResiduesByType(std::vector<Residue*> inputResidues,
                                                          cds::ResidueType queryType, const bool invert)
{ // Quality of life wrapper: calls the below function with one queryType.
    return selectResiduesByType(inputResidues, std::vector<cds::ResidueType> {queryType}, invert);
}

std::vector<Residue*> cdsSelections::selectResiduesByType(std::vector<Residue*> inputResidues,
                                                          std::vector<cds::ResidueType> queryTypes, const bool invert)
{
    std::vector<Residue*> selectedResidues;
    for (auto& residue : inputResidues)
    {
        bool contains = codeUtils::contains(queryTypes, residue->GetType());
        if ((contains && !invert) || (!contains && invert))
        {
            selectedResidues.push_back(residue);
        }
    }
    return selectedResidues;
}

// Not having Atom know which Residue it is in makes this funky. Make a decision about whether that happens or not.
Residue* cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(Residue* queryResidue,
                                                                    const std::string& queryAtomName)
{
    cds::Atom* queryAtom = queryResidue->FindAtom(queryAtomName);
    if (queryAtom == nullptr)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Warning: An atom named " + queryAtomName +
                      " not found in residue: " + cds::residueStringId(queryResidue));
        return nullptr;
    }
    cds::Atom* foreignAtomNeighbor = cdsSelections::selectNeighborNotInAtomVector(queryAtom, queryResidue->getAtoms());
    if (foreignAtomNeighbor == nullptr)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Warning: Did not find foreign neighbors for an atom named " + queryAtomName +
                      " in residue: " + cds::residueStringId(queryResidue));
        return nullptr;
    }
    for (auto& residueNeighbor : queryResidue->getNeighbors())
    {
        if (residueNeighbor->contains(foreignAtomNeighbor))
        {
            return residueNeighbor; // happy path.
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR,
              "Warning: Did not find a neighbor residue connected via " + queryAtomName +
                  " to residue: " + cds::residueStringId(queryResidue));
    return nullptr;
}

void cdsSelections::FindConnectedResidues(std::vector<Residue*>& visitedList, Residue* current)
{
    visitedList.push_back(current);
    for (auto& neighbor : current->getNeighbors())
    {
        if (!codeUtils::contains(visitedList, neighbor))
        {                                                                // Keep looking if neighbor wasn't yet visited.
            cdsSelections::FindConnectedResidues(visitedList, neighbor); // recursive function call
        }
    }
    return;
}

std::vector<Residue*> cdsSelections::selectResiduesWithinDistanceN(std::vector<Residue*> inputResidues,
                                                                   Residue* queryResidue, double queryDistance)
{
    std::vector<Residue*> foundResidues;
    cds::Coordinate queryCenter = cds::coordinateMean(cds::atomCoordinates(queryResidue->getAtoms()));
    for (auto& inputRes : inputResidues)
    {
        if (withinDistance(queryDistance, queryCenter, cds::coordinateMean(cds::atomCoordinates(inputRes->getAtoms()))))
        {
            foundResidues.push_back(inputRes);
        }
    }
    return foundResidues;
}
