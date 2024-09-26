#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <array>
#include <vector>

namespace cds
{
    template<typename T> void serializeNumbers(std::vector<T*> elements)
    {
        for (size_t n = 0; n < elements.size(); n++)
        {
            elements[n]->setNumber(n + 1);
        }
        return;
    }

    std::vector<int> atomNumbers(const std::vector<Atom*>& atoms);
    std::vector<std::string> atomNames(const std::vector<Atom*>& atoms);
    std::vector<std::string> atomElements(const std::vector<Atom*>& atoms);
    std::vector<int> atomAtomicNumbers(const std::vector<Atom*>& atoms);
    std::vector<std::string> atomTypes(const std::vector<Atom*>& atoms);
    std::vector<double> atomCharges(const std::vector<Atom*>& atoms);
    std::vector<int> residueNumbers(const std::vector<Residue*>& residues);
    std::vector<std::string> residueNames(const std::vector<Residue*>& residues);
    std::vector<ResidueType> residueTypes(const std::vector<Residue*>& residues);
    std::array<Atom*, 2> bondedAtomPair(std::vector<Atom*>& atomsA, std::vector<Atom*>& atomsB);
    std::vector<bool> atomsBondedTo(Atom* origin, std::vector<Atom*>& atoms);
} // namespace cds
#endif
