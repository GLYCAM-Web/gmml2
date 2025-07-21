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

    std::vector<Atom*> findCycleAtoms(Atom* const starterAtom)
    {
        glygraph::Graph<Atom> atomGraph(starterAtom);
        std::vector<std::pair<std::unordered_set<glygraph::Node<Atom>*>, std::unordered_set<glygraph::Edge<Atom>*>>>
            g1Cycles = cycle_decomp::totalCycleDetect(atomGraph);
        std::vector<Atom*> cycleAtoms;
        for (std::pair<std::unordered_set<glygraph::Node<Atom>*>, std::unordered_set<glygraph::Edge<Atom>*>>
                 currCyclePair : g1Cycles)
        {
            for (glygraph::Node<Atom>* currAtom : currCyclePair.first)
            {
                cycleAtoms.push_back(currAtom->getDerivedClass());
            }
        }
        return cycleAtoms;
    }

    // For now it's the lowest numbered (e.g. C1 lower than C6) ring atom connected to the ring oxygen.
    Atom* guessAnomericAtomByInternalNeighbors(const std::vector<Atom*> atoms)
    {
        if (atoms.empty())
        {
            throw std::runtime_error("Empty atom vector passed to guessAnomericAtomByInternalNeighbors");
        }
        std::vector<Atom*> cycleAtoms = findCycleAtoms(atoms.at(0));
        Atom* cycleOxygen = nullptr;
        for (auto& cycleAtom : cycleAtoms)
        {
            if (cycleAtom->getDerivedClass()->getElement() == "O")
            {
                cycleOxygen = cycleAtom->getDerivedClass();
            }
        }
        if (cycleOxygen == nullptr)
        {
            throw std::runtime_error(
                "Did not find a ring oxygen when trying to guess what the anomeric carbon is. Your "
                "atom names, are likely, strange.");
        }
        std::vector<Atom*> anomerCandidates = cycleOxygen->getNeighbors();
        // So deciding this isn't straightforward. For me use case of prep files this is fine. For reading form the PDB
        // I want a meeting first to figure out what we really need. Using a lambda to sort the candidates so that C1
        // appears before C2 etc.
        std::sort(
            anomerCandidates.begin(),
            anomerCandidates.end(),
            [](Atom*& a, Atom*& b) { return (a->getNumberFromName() < b->getNumberFromName()); });
        Atom* bestOne = anomerCandidates.front();
        return bestOne;
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

    std::vector<size_t> FindHeavyAtoms(const std::vector<Element>& elements)
    {
        std::vector<size_t> foundAtoms;
        foundAtoms.reserve(elements.size());
        std::vector<std::string> heavyList = {"C", "O", "N", "S", "P"};
        for (size_t n = 0; n < elements.size(); n++)
        {
            if (isHeavyElement(elements[n]))
            {
                foundAtoms.push_back(n);
            }
        }
        return foundAtoms;
    }
} // namespace gmml
