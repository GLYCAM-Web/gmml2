#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomCoordinates.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CodeUtils/logging.hpp"

using cds::Atom;
using cds::Molecule;
using cds::Residue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Molecule::Molecule(const std::string chainID) : Node<Molecule>(glygraph::invalid, {}), number_(0), chainId_(chainID)
{}

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

void Molecule::setResidues(std::vector<std::unique_ptr<Residue>> myResidues)
{ // This is good: myResidue contains a vector of unique_ptr, so you don't want to copy that.
    residues_ = std::move(myResidues);
}

Residue* Molecule::insertNewResidue(std::unique_ptr<Residue> myResidue, const Residue& positionReferenceResidue)
{
    auto position = this->findPositionOfResidue(&positionReferenceResidue);
    if (position != residues_.end())
    {
        ++position; // it is ok to insert at end(). I checked. It was ok. Ok.
        position = residues_.insert(position, std::move(myResidue));
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Could not create residue named " + myResidue->getName() + " as referenceResidue was not found\n");
    }
    return (*position).get(); // Dereference the reference to a uniquePtr, then use get() to create a raw ptr...
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

void Molecule::renumberResidues(int newStartNumber)
{
    for (auto& residue : this->getResidues())
    {
        residue->setNumber(newStartNumber++);
    }
}

//////////////////////////////////////////////////////////
//                    DISPLAY                           //
//////////////////////////////////////////////////////////
void Molecule::WritePdb(std::ostream& stream) const
{
    std::vector<Residue*> residues = getResidues();
    std::vector<ResidueType> types = residueTypes(residues);
    std::vector<bool> ter          = residueTER(types);
    PdbWriterData data             = toPdbWriterData(residues);
    cds::writeMoleculeToPdb(stream, codeUtils::indexVector(residues), ter, data);
    using cds::ResidueType; // to help readability of the Sugar, etc below
    std::vector<Residue*> selectedResidues =
        cdsSelections::selectResiduesByType(residues, {Sugar, Derivative, Aglycone, Undefined});
    std::vector<std::pair<Atom*, Atom*>> connectedAtoms = atomPairsConnectedToOtherResidues(selectedResidues);
    std::vector<std::pair<int, int>> connectionNumbers  = atomPairNumbers(connectedAtoms);
    cds::writeConectCards(stream, connectionNumbers);
}

void Molecule::WriteOff(std::ostream& stream)
{
    std::vector<Residue*> residues = getResidues();
    std::vector<Atom*> atoms       = getAtoms();
    cds::serializeNumbers(atoms);
    cds::serializeNumbers(residues);
    cds::WriteToOffFile(residues, stream, getName());
}
