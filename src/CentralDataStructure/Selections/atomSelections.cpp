#include "include/CentralDataStructure/Selections/atomSelections.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/templateGraph/TotalCycleDecomposition.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

namespace gmml
{
    Atom* getNonCarbonHeavyAtomNumbered(std::vector<Atom*> atoms, const std::string& queryNumber)
    {
        for (auto& atom : atoms)
        { // Assumes atom names like C2, O2 or N2. Nothing else should match.
            if (atom->getName().size() > 1)
            {
                const std::string number = atom->getName().substr(1);
                const std::string element = atom->getName().substr(0, 1); // Only the first character
                if ((number == queryNumber) && (element != "C") && (element != "H"))
                {
                    return atom;
                }
            }
        }
        return nullptr;
    }

    void FindConnectedAtoms(std::vector<Atom*>& visitedAtoms, Atom* currentAtom)
    {
        visitedAtoms.push_back(currentAtom);
        for (auto& neighbor : currentAtom->getNeighbors())
        {
            if (!util::contains(visitedAtoms, neighbor))
            {                                               // Keep looking if neighbor wasn't yet visited.
                FindConnectedAtoms(visitedAtoms, neighbor); // recursive function call
            }
        }
        return;
    }

    Atom* selectNeighborNotInAtomVector(const Atom* atomWithNeighbors, std::vector<Atom*> queryAtoms)
    {
        for (auto& neighbor : atomWithNeighbors->getNeighbors())
        {
            if (!util::contains(queryAtoms, neighbor))
            {
                return neighbor;
            }
        }
        return nullptr;
    }
} // namespace gmml
