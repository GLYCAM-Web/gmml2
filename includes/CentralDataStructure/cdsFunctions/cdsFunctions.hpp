#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP

#include "includes/CentralDataStructure/atom.hpp"

#include <array>
#include <vector>

namespace cds
{
    template<typename T> void serializeNumbers(std::vector<T*> elements)
    {
        unsigned int i = 0;
        for (auto& element : elements)
        {
            element->setNumber(++i);
        }
        return;
    }

    std::array<Atom*, 2> bondedAtomPair(std::vector<Atom*>& atomsA, std::vector<Atom*>& atomsB);
    std::vector<bool> atomsBondedTo(Atom* origin, std::vector<Atom*>& atoms);
} // namespace cds
#endif
