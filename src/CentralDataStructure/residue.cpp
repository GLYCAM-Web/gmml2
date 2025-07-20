#include "includes/CentralDataStructure/residue.hpp"

#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CodeUtils/biology.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/constants.hpp" // sNotSet
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <functional>
#include <string>
#include <vector>

using cds::Atom;
using cds::Residue;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
Residue::Residue(const std::string& residueName, const Residue* referenceResidue) : Node<Residue>(glygraph::invalid, {})
{
    this->setName(residueName);
    this->setNumber(referenceResidue->getNumber() + 1);
    this->determineType(residueName);
}

// Move Ctor.
Residue::Residue(Residue&& other) noexcept : Node<Residue>(other) { swap(*this, other); }

// Copy Ctor. Using copy-swap idiom. Call the base class copy ctor.
Residue::Residue(const Residue& other)
    : glygraph::Node<cds::Residue>(other), name_(other.name_), type_(other.type_), number_(other.number_)
{
    std::vector<cds::Atom*> otherAtoms = other.getAtoms();
    for (auto& atom : otherAtoms)
    { // create all the copies
        atoms_.push_back(std::make_unique<Atom>(*atom));
    }
    std::vector<cds::Atom*>::iterator ito = otherAtoms.begin();
    for (long unsigned int i = 0; i < otherAtoms.size(); i++)
    { // copy the connectivities
        std::vector<cds::Atom*> neighbors = otherAtoms.at(i)->getNeighbors();
        for (auto& neighbor : neighbors)
        {
            auto neighborPosition = std::find(otherAtoms.begin(), otherAtoms.end(), neighbor);
            if (neighborPosition != std::end(otherAtoms))
            { // neighbor could be in other residue, ignore here.
                auto difference = std::distance(ito, neighborPosition);
                int j = i + difference;
                std::string edgeName = atoms_.at(i)->getName() + "-" + atoms_.at(j)->getName();
                atoms_.at(i)->addNeighbor(edgeName, atoms_.at(j).get());
            }
        }
        ++ito;
    }
}

Residue& Residue::operator=(Residue other)
{
    this->Node<Residue>::operator=(other);
    swap(*this, other);
    return *this;
}

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////

std::vector<Atom*> Residue::getAtomsConditional(std::function<bool(Atom*)> condition) const
{
    std::vector<Atom*> atoms;
    atoms.reserve(atoms_.size());
    for (auto& atomPtr : atoms_)
    {
        Atom* atom = atomPtr.get();
        if (condition(atom))
        {
            atoms.push_back(atom);
        }
    }
    return atoms;
}

std::vector<Atom*> Residue::getAtoms() const
{
    auto all = [](Atom*) { return true; };
    return getAtomsConditional(all);
}

std::vector<Atom*> Residue::getHydrogenAtoms() const
{
    auto hydrogen = [](Atom* atom) { return atom->cachedElement() == MolecularMetadata::Element::H; };
    return getAtomsConditional(hydrogen);
}

std::vector<Atom*> Residue::getNonHydrogenAtoms() const
{
    auto nonHydrogen = [](Atom* atom) { return atom->cachedElement() != MolecularMetadata::Element::H; };
    return getAtomsConditional(nonHydrogen);
}

std::vector<Atom*> Residue::mutableAtoms()
{
    std::vector<Atom*> atoms;
    atoms.reserve(atoms_.size());
    for (auto& atomPtr : atoms_)
    {
        atoms.push_back(atomPtr.get());
    }
    return atoms;
}

const std::string Residue::GetParmName() const // If terminal, need to look up e.g. NPRO or CPRO instead of PRO.
{
    if (isNTerminal)
    {
        return "N" + this->getName();
    }
    else if (isCTerminal)
    {
        return "C" + this->getName();
    }
    return this->getName();
}

std::vector<std::string> Residue::getAtomNames() const
{
    std::vector<std::string> foundAtomNames;
    foundAtomNames.reserve(atoms_.size());
    for (auto& atom : atoms_)
    {
        foundAtomNames.push_back(atom.get()->getName());
    }
    return foundAtomNames;
}

//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
Atom* Residue::addAtom(std::unique_ptr<Atom> myAtom)
{
    atoms_.push_back(std::move(myAtom));
    return atoms_.back().get();
}

Atom* Residue::addAtomToFront(std::unique_ptr<Atom> myAtom)
{ // Yes expensive, but sometimes necessary.
    atoms_.insert(atoms_.begin(), std::move(myAtom));
    return atoms_.front().get();
}

bool Residue::moveAtomToLastPosition(const Atom* atom)
{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
    auto i = this->FindPositionOfAtom(atom); // auto makes my life easier
    if (i == atoms_.end())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Could not find atom in Residue to move to last position");
        return false; // atom not found maybe throw is better?
    }
    if (i == atoms_.end() - 1)
    {
        return true; // already in last position
    }
    std::iter_swap(i, atoms_.end() - 1);
    return true;
}

bool Residue::deleteAtom(const Atom* atom)
{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
    auto i = this->FindPositionOfAtom(atom); // auto makes my life easier
    if (i != atoms_.end())
    {
        // gmml::log(__LINE__, __FILE__, gmml::INF, "Atom " + atom->getName() + " has been erased. You're welcome.");
        i = atoms_.erase(i); // this costs a lot as everything after i gets shifted.
        return true;
    }
    return false;
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
typename std::vector<std::unique_ptr<Atom>>::iterator Residue::FindPositionOfAtom(const Atom* queryAtom)
{
    typename std::vector<std::unique_ptr<Atom>>::iterator i = atoms_.begin();
    typename std::vector<std::unique_ptr<Atom>>::iterator e = atoms_.end();
    while (i != e)
    {
        if (queryAtom == i->get())
        {
            return i;
        }
        else
        {
            ++i;
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Did not find " + queryAtom->getName() + " in atom records\n");
    return e;
}

Atom* Residue::FindAtom(const std::string queryName) const
{
    std::vector<Atom*> atoms = getAtoms();
    std::vector<std::string> names = atomNames(atoms);
    size_t index = codeUtils::indexOf(names, queryName);
    return index < atoms.size() ? atoms[index] : nullptr;
}

Atom* Residue::FindAtom(uint queryNumber) const
{
    std::vector<Atom*> atoms = getAtoms();
    std::vector<uint> numbers = atomNumbers(atoms);
    size_t index = codeUtils::indexOf(numbers, queryNumber);
    return index < atoms.size() ? atoms[index] : nullptr;
}

bool Residue::contains(const Atom* queryAtom) const
{
    std::vector<Atom*> atoms = this->getAtoms();
    return (std::find(atoms.begin(), atoms.end(), queryAtom) != atoms.end());
}

cds::ResidueType Residue::determineType(const std::string& residueName)
{
    if (codeUtils::contains(biology::proteinResidueNames, residueName))
    {
        this->SetType(ResidueType::Protein);
        return ResidueType::Protein;
    }
    // ToDo we want to figure out solvent, aglycone etc here too?.
    return ResidueType::Undefined;
}

std::string cds::residueStringId(cds::Residue* residue)
{
    std::function<std::string(const std::string&)> strOrNotSet = [](const std::string& str)
    { return str.empty() ? constants::sNotSet : str; };
    std::string name = residue->getName();
    std::string number = std::to_string(residue->getNumber());
    return codeUtils::join("_", codeUtils::vectorMap(strOrNotSet, {name, number, "", ""}));
}
