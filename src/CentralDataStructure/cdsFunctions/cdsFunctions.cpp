#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <array>
#include <vector>

std::vector<int> cds::atomNumbers(const std::vector<Atom*>& atoms)
{
    std::vector<int> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getNumber());
    }
    return result;
}

std::vector<std::string> cds::atomNames(const std::vector<Atom*>& atoms)
{
    std::vector<std::string> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getName());
    }
    return result;
}

std::vector<std::string> cds::atomElements(const std::vector<Atom*>& atoms)
{
    std::vector<std::string> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getElement());
    }
    return result;
}

std::vector<int> cds::atomAtomicNumbers(const std::vector<Atom*>& atoms)
{
    std::vector<int> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getAtomicNumber());
    }
    return result;
}

std::vector<std::string> cds::atomTypes(const std::vector<Atom*>& atoms)
{
    std::vector<std::string> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getType());
    }
    return result;
}

std::vector<double> cds::atomCharges(const std::vector<Atom*>& atoms)
{
    std::vector<double> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getCharge());
    }
    return result;
}

std::vector<int> cds::residueNumbers(const std::vector<Residue*>& residues)
{
    std::vector<int> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residue->getNumber());
    }
    return result;
}

std::vector<std::string> cds::residueNames(const std::vector<Residue*>& residues)
{
    std::vector<std::string> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residue->getName());
    }
    return result;
}

std::vector<cds::ResidueType> cds::residueTypes(const std::vector<Residue*>& residues)
{
    std::vector<ResidueType> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residue->GetType());
    }
    return result;
}

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
