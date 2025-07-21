#ifndef INCLUDES_READERS_PREP_PREPATOM_HPP
#define INCLUDES_READERS_PREP_PREPATOM_HPP

#include "include/CentralDataStructure/atom.hpp"
#include "include/util/constants.hpp"

#include <istream>
#include <ostream>
#include <string>

namespace gmml
{
    namespace prep
    { // repeated from common or utils as they should be here, not in gmml scope. Left over there as the old class needs
      // them until I delete it.

        enum TopologicalType
        {
            kTopTypeE,
            kTopTypeS,
            kTopTypeB,
            kTopType3,
            kTopType4,
            kTopTypeM
        };

        struct PrepAtomProperties
        {
            /*!< Sample line of the atom section of a prep file: 4 H1   H1  M  3  2  1  1.000    90.0     180.0     0.0
             */
            std::string type = ""; /*!< Atom type; fill by the third column of the residue section of the file */
            TopologicalType topologicalType = kTopTypeM; /*!< Topological type (for chain extraction of the residue);
                                                            fill by th 4th column of the residue section of the file */
            int bondIndex = 0;     /*!< Bond index; fill by the 5th column of the residue section of the file */
            int angleIndex = 0;    /*!< Angle index; fill by the 6th column of the residue section of the file */
            int dihedralIndex = 0; /*!< Dihedral index; fill by the 7th column of the residue section of the file */
            double bondLength =
                constants::dNotSet; /*!< Bond length; fill by the 8th column of the residue section of the file */
            double angle = constants::dNotSet; /*!< Angle; fill by the 9th column of the residue section of the file */
            double dihedral =
                constants::dNotSet; /*!< Dihedral; fill by the 10th column of the residue section of the file */
            int visitCount = 0;
        };

        class PrepAtom : public Atom
        {
          public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PrepAtom(const std::string& line);
            PrepAtom() = default;
            ~PrepAtom() = default;
            PrepAtom(PrepAtom&& other) noexcept;               // Move Ctor
            explicit PrepAtom(const PrepAtom& other) noexcept; // Copy Ctor
            PrepAtom& operator=(PrepAtom other) noexcept;      // Move and Copy assignment operator

            friend void swap(PrepAtom& lhs, PrepAtom& rhs) noexcept
            {
                using std::swap;
                swap(static_cast<Atom&>(lhs), static_cast<Atom&>(rhs));
                swap(lhs.properties, rhs.properties);
            }

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            void Determine3dCoordinate();

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out) const;
            void Write(std::ostream& stream) const;

            PrepAtomProperties properties;
        };
    } // namespace prep
} // namespace gmml

#endif
