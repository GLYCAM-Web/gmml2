#ifndef INCLUDE_READERS_PREP_PREPATOM_HPP
#define INCLUDE_READERS_PREP_PREPATOM_HPP

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

        static const std::vector<std::string> topologicalTypeNames {"E", "S", "B", "3", "4", "M"};

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
        };

        void initializePrepAtom(Atom* atom, PrepAtomProperties& properties, const std::string& line);
        void determine3dCoordinate(Atom* atom, const PrepAtomProperties& properties);
        void print(Atom* atom, const PrepAtomProperties& properties, std::ostream& out);
        void write(Atom* atom, const PrepAtomProperties& properties, std::ostream& stream);
    } // namespace prep
} // namespace gmml

#endif
