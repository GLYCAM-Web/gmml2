#include "includes/CentralDataStructure/Selections/atomSelections.hpp"

#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/MolecularModeling/TemplateGraph/Algorithms/include/TotalCycleDecomposition.hpp"

Atom* cdsSelections::getNonCarbonHeavyAtomNumbered(std::vector<Atom*> atoms, const std::string& queryNumber)
{
    for (auto& atom : atoms)
    { // Assumes atom names like C2, O2 or N2. Nothing else should match.
        if (atom->getName().size() > 1)
        {
            const std::string number  = atom->getName().substr(1);
            const std::string element = atom->getName().substr(0, 1); // Only the first character
            if ((number == queryNumber) && (element != "C") && (element != "H"))
            {
                return atom;
            }
        }
    }
    return nullptr;
}

void cdsSelections::FindConnectedAtoms(std::vector<Atom*>& visitedAtoms, Atom* currentAtom)
{
    visitedAtoms.push_back(currentAtom);
    for (auto& neighbor : currentAtom->getNeighbors())
    {
        if (!codeUtils::contains(visitedAtoms, neighbor))
        {                                                              // Keep looking if neighbor wasn't yet visited.
            cdsSelections::FindConnectedAtoms(visitedAtoms, neighbor); // recursive function call
        }
    }
    return;
}

std::vector<Atom*> cdsSelections::findCycleAtoms(cds::Atom* const starterAtom)
{
    glygraph::Graph<cds::Atom> atomGraph(starterAtom);
    std::vector<std::pair<std::unordered_set<glygraph::Node<Atom>*>, std::unordered_set<glygraph::Edge<Atom>*>>>
        g1Cycles = cycle_decomp::totalCycleDetect(atomGraph);
    std::vector<Atom*> cycleAtoms;
    for (std::pair<std::unordered_set<glygraph::Node<Atom>*>, std::unordered_set<glygraph::Edge<Atom>*>> currCyclePair :
         g1Cycles)
    {
        for (glygraph::Node<Atom>* currAtom : currCyclePair.first)
        {
            cycleAtoms.push_back(currAtom->getDerivedClass());
        }
    }
    return cycleAtoms;
}

// For now it's the lowest numbered (e.g. C1 lower than C6) ring atom connected to the ring oxygen.
Atom* cdsSelections::guessAnomericAtomByInternalNeighbors(const std::vector<cds::Atom*> atoms)
{
    if (atoms.empty())
    {
        throw std::runtime_error("Empty atom vector passed to guessAnomericAtomByInternalNeighbors");
    }
    std::vector<Atom*> cycleAtoms = cdsSelections::findCycleAtoms(atoms.at(0));
    cds::Atom* cycleOxygen        = nullptr;
    for (auto& cycleAtom : cycleAtoms)
    {
        if (cycleAtom->getDerivedClass()->getElement() == "O")
        {
            cycleOxygen = cycleAtom->getDerivedClass();
        }
    }
    if (cycleOxygen == nullptr)
    {
        throw std::runtime_error("Did not find a ring oxygen when trying to guess what the anomeric carbon is. Your "
                                 "atom names, are likely, strange.");
    }
    std::vector<Atom*> anomerCandidates = cycleOxygen->getNeighbors();
    // So deciding this isn't straightforward. For me use case of prep files this is fine. For reading form the PDB I
    // want a meeting first to figure out what we really need. Using a lambda to sort the candidates so that C1 appears
    // before C2 etc.
    std::sort(anomerCandidates.begin(), anomerCandidates.end(),
              [](Atom*& a, Atom*& b)
              {
                  return (a->getNumberFromName() < b->getNumberFromName());
              });
    Atom* bestOne = anomerCandidates.front();
    return bestOne;
}

Atom* cdsSelections::selectNeighborNotInAtomVector(const Atom* atomWithNeighbors, std::vector<Atom*> queryAtoms)
{
    for (auto& neighbor : atomWithNeighbors->getNeighbors())
    {
        if (!codeUtils::contains(queryAtoms, neighbor))
        {
            return neighbor;
        }
    }
    return nullptr;
}

std::vector<size_t> cdsSelections::FindHeavyAtoms(const std::vector<MolecularMetadata::Element>& elements)
{
    std::vector<size_t> foundAtoms;
    foundAtoms.reserve(elements.size());
    std::vector<std::string> heavyList = {"C", "O", "N", "S", "P"};
    for (size_t n = 0; n < elements.size(); n++)
    {
        if (MolecularMetadata::isHeavyElement(elements[n]))
        {
            foundAtoms.push_back(n);
        }
    }
    return foundAtoms;
}

unsigned long int cdsSelections::CountAtomsWithinBondingDistance(const Atom* queryAtom, std::vector<Atom*> otherAtoms)
{
    unsigned long int count = 0;
    for (auto& otherAtom : otherAtoms)
    {
        if (isWithinBondingDistance(queryAtom, otherAtom))
        {
            ++count;
        }
    }
    return count;
}
