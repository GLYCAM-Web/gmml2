#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

using cds::Atom;
using cds::Molecule;
using cds::Residue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Molecule::Molecule(std::vector<Residue*>& residues) : Node<Molecule>(glygraph::invalid, {}), number_(0)
{
    for (auto& residue : residues)
    {
        residues_.push_back(std::make_unique<Residue>(std::move(*residue)));
    }
}

// Move Ctor
Molecule::Molecule(Molecule&& other) noexcept : glygraph::Node<cds::Molecule>(other)
{
    residues_ = std::move(other.residues_);
    number_   = std::move(other.number_);
}

// Copy Ctor
Molecule::Molecule(const Molecule& other) : glygraph::Node<cds::Molecule>(other), number_(other.number_)
{
    for (auto& residue : other.residues_)
    {
        residues_.push_back(std::make_unique<Residue>((*residue.get())));
    }
}

// Move and Copy assignment operator
Molecule& Molecule::operator=(Molecule other)
{
    glygraph::Node<cds::Molecule>::operator=(other); // ToDo ok?
    swap(*this, other);
    return *this;
}

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<Residue*> Molecule::getResidues() const
{
    std::vector<Residue*> residues;
    residues.reserve(residues_.size());
    for (auto& residuePtr : residues_)
    {
        residues.push_back(residuePtr.get());
    }
    return residues;
}

std::vector<Atom*> Molecule::getAtoms() const
{
    std::vector<Atom*> atoms;
    for (auto& residue : residues_)
    {
        codeUtils::insertInto(atoms, residue->getAtoms());
    }
    return atoms;
}

std::vector<Atom*> Molecule::mutableAtoms()
{
    std::vector<Atom*> atoms;
    for (auto& residue : residues_)
    {
        codeUtils::insertInto(atoms, residue->mutableAtoms());
    }
    return atoms;
}

//////////////////////////////////////////////////////////
//                    MUTATORS                          //
//////////////////////////////////////////////////////////
void Molecule::swapResiduePosition(Residue* queryResidue, size_t newPosition)
{
    for (size_t oldPosition = 0; oldPosition < residues_.size(); oldPosition++)
    {
        if (residues_[oldPosition].get() == queryResidue)
        {
            if (oldPosition != newPosition)
            {
                std::swap(residues_[oldPosition], residues_[newPosition]);
            }
            break;
        }
    }
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
Residue* Molecule::addResidue(std::unique_ptr<Residue> myResidue)
{ // This is good: myResidue contains a vector of unique_ptr, so you don't want to copy that.
    residues_.push_back(std::move(myResidue));
    return residues_.back().get();
}

Residue* Molecule::insertNewResidue(std::unique_ptr<Residue> myResidue, size_t position)
{
    if (position <= residues_.size())
    {
        auto it = residues_.insert(residues_.begin() + position, std::move(myResidue));
        return it->get();
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Could not create residue named " + myResidue->getName() + " as referenceResidue was not found\n");
    }
    return nullptr;
}

Residue* Molecule::insertNewResidue(std::unique_ptr<Residue> myResidue, const Residue& positionReferenceResidue)
{
    auto position = this->findPositionOfResidue(&positionReferenceResidue);
    return insertNewResidue(std::move(myResidue), 1 + position - residues_.begin());
}

std::vector<std::unique_ptr<Residue>>::iterator Molecule::findPositionOfResidue(const Residue* queryResidue)
{
    auto it = std::find_if(residues_.begin(), residues_.end(),
                           [&](auto& i)
                           {
                               return queryResidue == i.get();
                           });
    if (it == residues_.end())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Did not find position of " + queryResidue->getName() +
                      " in vector\n"); // every class should have a print?
    }
    return it;
}

std::vector<Residue*> Molecule::getResidues(std::vector<std::string> queryNames)
{
    return codeUtils::getElementsWithNames(this->getResidues(), queryNames);
}

Residue* Molecule::getResidue(const std::string& queryName)
{
    return codeUtils::findElementWithName(this->getResidues(), queryName);
}

void Molecule::deleteResidue(Residue* residue)
{
    auto i = this->findPositionOfResidue(residue); // auto makes my life easier
    if (i != residues_.end())
    {
        i = residues_.erase(i);
    }
    return;
}
