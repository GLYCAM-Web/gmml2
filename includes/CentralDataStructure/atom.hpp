#ifndef INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CodeUtils/constants.hpp" // dNotSet
#include "includes/CodeUtils/references.hpp"

#include <string>
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
        }

        //////////////////////////////////////////////////////////
        //                       ACCESSORS                      //
        //////////////////////////////////////////////////////////
        Coordinate coordinate() const;
        CoordinateReference coordinateReference();
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

        inline bool isVisible() const
        {
            return isVisible_;
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

        inline void makeInvisible()
        {
            isVisible_ = false;
        }

        void setCoordinate(const Coordinate& c);
        void setCurrentCoordinate(size_t coordinateIndex);
        void addCoordinate(const Coordinate& c);
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

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        size_t currentCoordinate_ = -1;
        std::vector<Coordinate> allCoordinates_;
        double charge_                      = constants::dNotSet;
        std::string atomType_               = " ";
        int number_                         = constants::iNotSet;
        MolecularMetadata::Element element_ = MolecularMetadata::Element::Unknown;
        bool gotElement_                    = false;
        bool isVisible_                     = true;
    };
} // namespace cds
#endif
