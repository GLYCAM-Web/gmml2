#include "include/CentralDataStructure/atom.hpp"

#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/logging.hpp"

#include <ctype.h> // isalpha
#include <sstream>

namespace gmml
{
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    //////////////////////////////////////////////////////////
    Atom::Atom(const std::string name, const Coordinate& coord) : Node<Atom>(name, {})
    {
        this->setCoordinate(coord);
        this->setNumber(1); // Seems like a fine default?
    }

    // Move Ctor
    Atom::Atom(Atom&& other) noexcept : Node<Atom>(other) { swap(*this, other); }

    // Copy Ctor
    Atom::Atom(const Atom& other) noexcept
        : glygraph::Node<Atom>(other), coordinate_(other.coordinate_), charge_(other.charge_),
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
        const std::string name = this->getName();
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
    void Atom::setElement(Element element)
    {
        element_ = element;
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
        util::log(__LINE__, __FILE__, util::WAR, "Did not find an element for atom named: " + name);
        return "";
    }

    Element Atom::cachedElement()
    {
        if (!gotElement_)
        {
            element_ = toElement(getElement());
            gotElement_ = true;
        }
        return element_;
    }

    int Atom::getAtomicNumber() const { return findElementAtomicNumber(this->getElement()); }

    std::string Atom::getId() const { return this->getName() + "_" + std::to_string(this->getIndex()); }

    //////////////////////////////////////////////////////////
    //                   OVERLOADED OPERATORS               //
    //////////////////////////////////////////////////////////
    bool Atom::operator==(const Atom& otherAtom) { return (this->getIndex() == otherAtom.getIndex()); }

    bool Atom::operator!=(const Atom& otherAtom) { return !(operator==(otherAtom)); }
} // namespace gmml
