#include "includes/CentralDataStructure/residue.hpp"

#include "../../includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "../../includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/constants.hpp" // sNotSet
#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
#include "includes/CodeUtils/biology.hpp"

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
Residue::Residue(Residue&& other) noexcept : Node<Residue>(other)
{
    swap(*this, other);
}

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
                auto difference      = std::distance(ito, neighborPosition);
                int j                = i + difference;
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

std::vector<Atom*> Residue::getAtoms() const
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
    if (this->containsLabel("NTerminal"))
    {
        return "N" + this->getName();
    }
    else if (this->containsLabel("CTerminal"))
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

std::string Residue::getStringId(std::string moleculeNumber) const
{
    const std::string insertionCode = constants::sNotSet;
    pdb::ResidueId temp(this->getName(), std::to_string(this->getNumber()), insertionCode, moleculeNumber);
    return temp.print();
}

std::vector<Coordinate*> Residue::getCoordinates() const
{
    return cds::getCoordinatesFromAtoms(this->getAtoms());
}

void Residue::insertCoordinatesInto(std::vector<Coordinate>& coordinates) const
{
    for (auto& a : atoms_)
    {
        coordinates.push_back(*a.get()->getCoordinate());
    }
}

pdb::ResidueId Residue::getId() const
{
    pdb::ResidueId temp(this->getName(), std::to_string(this->getNumber()), constants::sNotSet, constants::sNotSet);
    return temp;
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

// void Residue::addAtom(Atom* myAtom)
//{
//     atoms_.push_back(std::make_unique<Atom>(myAtom));
//     return;
// }

bool Residue::moveAtomToLastPosition(const Atom* atom)
{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
    auto i = this->FindPositionOfAtom(atom); // auto makes my life easier
    if (*i == atoms_.back())
    {
        return true; // already in last position
    }
    if (i == atoms_.end())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Could not find atom in Residue to move to last position");
        return false; // atom not found maybe throw is better?
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
    return codeUtils::findElementWithName(this->getAtoms(), queryName);
}

Atom* Residue::FindAtom(const int& queryNumber) const
{
    return codeUtils::findElementWithNumber(this->getAtoms(), queryNumber);
}

bool Residue::contains(const Atom* queryAtom) const
{
    std::vector<Atom*> atoms = this->getAtoms();
    return (std::find(atoms.begin(), atoms.end(), queryAtom) != atoms.end());
}

std::vector<const Atom*> Residue::getAtomsConnectedToOtherResidues() const
{
    std::vector<const Atom*> foundAtoms;
    std::vector<Atom*> residueAtoms = this->getAtoms();
    for (auto& atom : residueAtoms)
    {
        for (auto& neighbor : atom->getNeighbors())
        { // check if neighbor is not one of the atoms in this residue.
            if (std::find(residueAtoms.begin(), residueAtoms.end(), neighbor) == residueAtoms.end())
            {
                foundAtoms.push_back(atom);
            }
        }
    }
    return foundAtoms;
}

void Residue::findAtomPairsConnectedToOtherResidues(std::vector<std::pair<const Atom*, const Atom*>>& foundAtoms) const
{
    // std::vector<std::pair<const Atom*, const Atom*>> foundAtoms;
    std::vector<Atom*> residueAtoms = this->getAtoms();
    for (auto& atom : residueAtoms)
    { // only "child" neighbors or we find same pair twice
        for (auto& neighbor : atom->getChildren())
        { // check if neighbor is not one of the atoms in this residue.
            if (std::find(residueAtoms.begin(), residueAtoms.end(), neighbor) == residueAtoms.end())
            {
                foundAtoms.push_back({atom, neighbor});
            }
        }
    }
    return;
}

void Residue::MakeDeoxy(const std::string oxygenNumber)
{ // if oxygenNumber is 6, then C6-O6-H6O becomes C6-Hd
    Atom* hydrogenAtom = this->FindAtom("H" + oxygenNumber + "O");
    Atom* oxygenAtom   = this->FindAtom("O" + oxygenNumber);
    Atom* carbonAtom   = this->FindAtom("C" + oxygenNumber);
    // Add O and H charge to the C atom.
    carbonAtom->setCharge(carbonAtom->getCharge() + oxygenAtom->getCharge() + hydrogenAtom->getCharge());
    // Delete the H of O-H
    this->deleteAtom(hydrogenAtom);
    // Now transform the Oxygen to a Hd. Easier than deleting O and creating H. Note: this H looks weird in LiteMol as
    // bond length is too long.
    //    std::string newID = oxygenAtom->getId();
    //    newID.replace(0,oxygenAtom->getName().size(),"Hd");
    oxygenAtom->setName("Hd");
    oxygenAtom->setType("H1");
    oxygenAtom->setCharge(0.0000);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Completed MakeDeoxy\n");
}

cds::ResidueType Residue::determineType(const std::string& residueName)
{
    if (std::find(biology::proteinResidueNames.begin(), biology::proteinResidueNames.end(), residueName) !=
        biology::proteinResidueNames.end())
    {
        this->SetType(ResidueType::Protein);
        return ResidueType::Protein;
    }
    // ToDo we want to figure out solvent, aglycone etc here too?.
    return ResidueType::Undefined;
}

//////////////////////////////////////////////////////////
//                    DISPLAY                           //
//////////////////////////////////////////////////////////
void Residue::Print(std::ostream& out) const
{
    out << "Name: " << this->getName() << ", index:" << this->getIndex() << ":\n";
    for (auto& atom : this->getAtoms())
    {
        atom->Print(out);
        out << "\n";
    }
    out << "\n";
}
