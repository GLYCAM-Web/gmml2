#include "includes/CentralDataStructure/atom.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

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
    : glygraph::Node<cds::Atom>(other), coordinate_(other.coordinate_), charge_(other.charge_),
      atomType_(other.atomType_), number_(other.number_)
{}

// Move and Copy assignment operator
Atom& Atom::operator=(Atom other) noexcept
{
    this->Node<Atom>::operator=(other);
    swap(*this, other);
    return *this;
}

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////

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
void Atom::setElement(MolecularMetadata::Element element)
{
    element_    = element;
    gotElement_ = true;
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
    return !(operator==(otherAtom));
}
