#ifndef INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CodeUtils/constants.hpp" // dNotSet

#include <string>
#include <iostream>
#include <vector>
#include <memory> // unique_ptr

namespace cds
{
    class Atom : public glygraph::Node<Atom>
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTORS                   //
        //////////////////////////////////////////////////////////
        Atom() : Node<Atom>(glygraph::invalid, {}) {};
        Atom(const std::string name, const Coordinate& coord);
        Atom(Atom&& other) noexcept;               // Move Ctor
        explicit Atom(const Atom& other) noexcept; // Copy Ctor
        Atom& operator=(Atom other) noexcept;      // Move and Copy assignment operator

        friend void swap(Atom& lhs, Atom& rhs) noexcept
        {
            using std::swap;
            swap(lhs.currentCoordinate_, rhs.currentCoordinate_);
            swap(lhs.allCoordinates_, rhs.allCoordinates_);
            swap(lhs.charge_, rhs.charge_);
            swap(lhs.atomType_, rhs.atomType_);
            swap(lhs.number_, rhs.number_);
            std::cout << "Swap triggered for Atom, everything in " << rhs.getName() << " swapped with " << lhs.getName()
                      << std::endl;
        }

        //////////////////////////////////////////////////////////
        //                       ACCESSORS                      //
        //////////////////////////////////////////////////////////
        Coordinate coordinate() const;
        Coordinate* coordinatePointer();
        unsigned int getNumberOfCoordinateSets() const;

        inline double getCharge() const
        {
            return charge_;
        }

        inline std::string getType() const
        {
            return atomType_;
        }

        inline int getNumber() const
        {
            return number_;
        }

        unsigned int getNumberFromName() const;

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void setCharge(const double c)
        {
            charge_ = c;
        }

        inline void setType(const std::string s)
        {
            atomType_ = s;
        }

        inline void setNumber(const int i)
        {
            number_ = i;
        }

        void setCoordinate(const Coordinate& c);
        void setCurrentCoordinate(unsigned int coordinateIndex = 0);
        Coordinate* addCoordinate(const Coordinate& c);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string getElement() const;
        MolecularMetadata::Element cachedElement();
        int getAtomicNumber() const;
        virtual std::string getId() const;
        //////////////////////////////////////////////////////////
        //                   OVERLOADED OPERATORS               //
        //////////////////////////////////////////////////////////
        bool operator==(const Atom& otherAtom);
        bool operator!=(const Atom& otherAtom);
        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        virtual void Print(std::ostream& out) const;

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        Coordinate* currentCoordinate_ = nullptr;
        std::vector<std::unique_ptr<Coordinate>> allCoordinates_;
        double charge_                      = constants::dNotSet;
        std::string atomType_               = " ";
        int number_                         = constants::iNotSet;
        MolecularMetadata::Element element_ = MolecularMetadata::Element::Unknown;
        bool gotElement_                    = false;
    };
} // namespace cds
#endif
