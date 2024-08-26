#include "includes/CentralDataStructure/atom.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"
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

MolecularMetadata::Element Atom::cachedElement()
{
    if (!gotElement_)
    {
        element_    = MolecularMetadata::toElement(getElement());
        gotElement_ = true;
    }
    return element_;
}

int Atom::getAtomicNumber() const
{
    return MolecularMetadata::findElementAtomicNumber(this->getElement());
}

std::string Atom::getId() const
{
    return this->getName() + "_" + std::to_string(this->getIndex());
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

cds::Sphere cds::coordinateWithRadius(Atom* atom)
{
    auto element = atom->cachedElement();
    return {MolecularMetadata::vanDerWaalsRadius(element), *atom->getCoordinate()};
}

std::vector<Coordinate*> cds::getCoordinatesFromAtoms(std::vector<cds::Atom*> atoms)
{
    std::vector<Coordinate*> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        coordinates.push_back(atom->getCoordinate());
    }
    return coordinates;
}

std::vector<cds::Sphere> cds::getCoordinatesWithRadiiFromAtoms(std::vector<cds::Atom*> atoms)
{
    std::vector<cds::Sphere> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        coordinates.push_back(coordinateWithRadius(atom));
    }
    return coordinates;
}
