#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/atom.hpp"

#include <array>
#include <vector>

std::array<cds::Atom*, 2> cds::bondedAtomPair(std::vector<Atom*>& atomsA, std::vector<Atom*>& atomsB)
{
    for (auto& a : atomsA)
    {
        for (auto& b : atomsB)
        {
            if (a->isNeighbor(b))
            {
                return {a, b};
            }
        }
    }
    return {nullptr, nullptr};
}

std::vector<bool> cds::atomsBondedTo(Atom* origin, std::vector<Atom*>& atoms)
{
    std::vector<bool> result;
    for (size_t n = 0; n < atoms.size(); n++)
    {
        result.push_back((atoms[n] == origin) || origin->isNeighbor(atoms[n]));
    }
    return result;
}
