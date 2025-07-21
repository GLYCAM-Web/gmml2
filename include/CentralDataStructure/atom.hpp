#ifndef INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/templateGraph/Node.hpp"
#include "include/util/constants.hpp"

#include <memory> // unique_ptr
#include <string>
#include <vector>

namespace gmml
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
            swap(lhs.coordinate_, rhs.coordinate_);
            swap(lhs.charge_, rhs.charge_);
            swap(lhs.atomType_, rhs.atomType_);
            swap(lhs.number_, rhs.number_);
        }

        //////////////////////////////////////////////////////////
        //                       ACCESSORS                      //
        //////////////////////////////////////////////////////////
        inline Coordinate coordinate() const { return coordinate_; }

        inline double getCharge() const { return charge_; }

        inline std::string getType() const { return atomType_; }

        inline uint getNumber() const { return number_; }

        inline bool isVisible() const { return isVisible_; }

        unsigned int getNumberFromName() const;

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void setCharge(const double c) { charge_ = c; }

        inline void setType(const std::string s) { atomType_ = s; }

        inline void setNumber(uint i) { number_ = i; }

        inline void makeInvisible() { isVisible_ = false; }

        void setElement(Element element);

        inline void setCoordinate(const Coordinate& c) { coordinate_ = c; }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string getElement() const;
        Element cachedElement();
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
        Coordinate coordinate_ = {0.0, 0.0, 0.0};
        double charge_ = constants::dNotSet;
        std::string atomType_ = " ";
        uint number_ = constants::iNotSet;
        Element element_ = Element::Unknown;
        bool gotElement_ = false;
        bool isVisible_ = true;
    };
} // namespace gmml
#endif
