#include "include/CentralDataStructure/Selections/residueSelections.hpp"

#include "include/CentralDataStructure/Selections/atomSelections.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/readers/Pdb/pdbData.hpp"
#include "include/readers/Pdb/pdbResidue.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

namespace gmml
{
    std::vector<Residue*> selectResiduesByType(
        std::vector<Residue*> inputResidues, ResidueType queryType, const bool invert)
    { // Quality of life wrapper: calls the below function with one queryType.
        return selectResiduesByType(inputResidues, std::vector<ResidueType> {queryType}, invert);
    }

    std::vector<Residue*> selectResiduesByType(
        std::vector<Residue*> inputResidues, std::vector<ResidueType> queryTypes, const bool invert)
    {
        std::vector<Residue*> selectedResidues;
        for (auto& residue : inputResidues)
        {
            bool contains = util::contains(queryTypes, residue->GetType());
            if ((contains && !invert) || (!contains && invert))
            {
                selectedResidues.push_back(residue);
            }
        }
        return selectedResidues;
    }

    void FindConnectedResidues(std::vector<Residue*>& visitedList, Residue* current)
    {
        visitedList.push_back(current);
        for (auto& neighbor : current->getNeighbors())
        {
            if (!util::contains(visitedList, neighbor))
            {                                                 // Keep looking if neighbor wasn't yet visited.
                FindConnectedResidues(visitedList, neighbor); // recursive function call
            }
        }
        return;
    }
} // namespace gmml
