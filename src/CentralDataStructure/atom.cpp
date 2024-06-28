#include "includes/CentralDataStructure/atom.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp" // bondIfClose
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/MolecularMetadata/atomicNumbers.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"

#include <ctype.h> // isalpha
#include <sstream>

using cds::Atom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Atom::Atom(const std::string name, const Coordinate& coord) : Node<Atom>(name, {})
{
    this->setCoordinate(coord);
    this->setNumber(1); // Seems like a fine default?
}

// Move Ctor
Atom::Atom(Atom&& other) noexcept : Node<Atom>(other)
{
    swap(*this, other);
}

// Copy Ctor
Atom::Atom(const Atom& other) noexcept
    : glygraph::Node<cds::Atom>(other), currentCoordinate_(other.currentCoordinate_), charge_(other.charge_),
      atomType_(other.atomType_), number_(other.number_)
{
    //    gmml::log(__LINE__, __FILE__, gmml::INF,
    //              "Atom copy ctor creating " + this->getName() + "_" + std::to_string(this->getIndex()));
    for (auto& coord : other.allCoordinates_)
    {
        allCoordinates_.push_back(std::make_unique<Coordinate>((*coord.get())));
    }
    currentCoordinate_ = allCoordinates_.at(0).get();
}

// Move and Copy assignment operator
Atom& Atom::operator=(Atom other) noexcept
{
    this->Node<Atom>::operator=(other);
    swap(*this, other);
    currentCoordinate_ = allCoordinates_.at(0).get();
    return *this;
}

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
Coordinate* Atom::getCoordinate()
{
    return currentCoordinate_;
}

Coordinate* Atom::getCoordinate() const
{
    return currentCoordinate_;
}

unsigned int Atom::getNumberOfCoordinateSets() const
{
    return allCoordinates_.size();
}

unsigned int Atom::getNumberFromName() const
{
    const std::string name        = this->getName();
    std::string convertableNumber = "";
    for (size_t i = 1; i < name.size(); i++)
    {
        if (isdigit(name[i]))
        {
            convertableNumber += name[i];
        }
    }
    if (convertableNumber.empty())
    {
        return 0;
    }
    return std::stoi(convertableNumber);
}

//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
void Atom::setCoordinate(const Coordinate& newCoord)
{
    currentCoordinate_ = this->addCoordinate(newCoord);
}

void Atom::setCurrentCoordinate(unsigned int coordinateIndex)
{
    if (allCoordinates_.size() <= coordinateIndex)
    {
        std::stringstream ss;
        ss << "Error: requested coordinateIndex: " << coordinateIndex
           << " that doesn't exist in Atom class as allCoordinates_ size is " << allCoordinates_.size();
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    currentCoordinate_ = allCoordinates_.at(coordinateIndex).get();
}

Coordinate* Atom::addCoordinate(const Coordinate& newCoord)
{
    allCoordinates_.push_back(std::make_unique<Coordinate>(newCoord));
    return allCoordinates_.back().get();
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void Atom::addBond(Atom* otherAtom)
{
    this->addNeighbor("atomicBond", otherAtom);
    return;
}

void Atom::bondIfClose(Atom* otherAtom)
{
    if (this->isWithinBondingDistance(otherAtom))
    {
        this->addBond(otherAtom);
        // std::stringstream ss;
        // std::cout << "Bonded " << this->getName() << "_" << this->getIndex() << " to " << otherAtom->getName() << "_"
        // << otherAtom->getIndex() << "\n"; gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
    }
    return;
}

std::string Atom::getElement() const // derived classes should overwrite if more explicit about element.
{
    std::string name = this->getName();
    if (!name.empty())
    {
        if (isalpha(name.at(0))) // if first char is in the alphabet
        {
            return name.substr(0, 1); // return first character as string
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR, "Did not find an element for atom named: " + name);
    return "";
}

int Atom::getAtomicNumber() const
{
    return MolecularMetadata::findElementAtomicNumber(this->getElement());
}

std::string Atom::getId() const
{
    return this->getName() + "_" + std::to_string(this->getIndex());
}

bool Atom::isWithinBondingDistance(const Atom* otherAtom) const
{
    double maxLength = atomicBonds::getMaxBondLengthByAtomType(this->getElement(), otherAtom->getElement());
    return withinDistance(maxLength, *this->getCoordinate(), *otherAtom->getCoordinate());
}

double Atom::calculateDistance(const Atom* otherAtom) const
{
    return distance(*this->getCoordinate(), *otherAtom->getCoordinate());
}

//////////////////////////////////////////////////////////
//                   OVERLOADED OPERATORS               //
//////////////////////////////////////////////////////////
bool Atom::operator==(const Atom& otherAtom)
{
    return (this->getIndex() == otherAtom.getIndex());
}

bool Atom::operator!=(const Atom& otherAtom)
{
    return (this->getIndex() != otherAtom.getIndex());
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void Atom::Print(std::ostream& out) const
{
    out << this->getName() << ", ";
    return;
}
