#include "includes/CentralDataStructure/Selections/residueSelections.hpp"

#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
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
