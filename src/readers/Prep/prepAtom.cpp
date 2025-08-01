#include "include/readers/Prep/prepAtom.hpp"

#include "include/geometry/measurements.hpp" //get_cartesian_point_from_internal_coords()
#include "include/util/casting.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <iomanip>
#include <ios>
#include <sstream>
#include <stdexcept>

namespace gmml
{
    namespace prep
    {
        namespace
        {
            TopologicalType extractAtomTopologicalType(std::istream& ss)
            {
                std::string s;
                ss >> s;
                size_t index = util::indexOf(topologicalTypeNames, s);
                if (index == topologicalTypeNames.size())
                {
                    throw std::runtime_error("Unknown topological type");
                }
                return TopologicalType(index);
            }

            std::vector<Atom*> findDihedralAtoms(Atom* initial)
            {
                std::vector<Atom*> result = {initial};
                for (int currentDepth = 0; currentDepth < 3; currentDepth++)
                {
                    Atom* parent = result.back()->getParents().front(); // Go up the first parent only. Loops may create
                                                                        // another parent, but they should be ignored.
                    result.push_back(parent);
                }
                return result;
            }
        } // namespace

        void initializePrepAtom(Atom* atom, PrepAtomProperties& properties, const std::string& line)
        {
            std::stringstream ss(line);
            atom->setNumber(util::extractFromStream(ss, int()));
            atom->setName(util::extractFromStream(ss, std::string()));
            atom->setType(util::extractFromStream(ss, std::string()));
            properties.topologicalType = extractAtomTopologicalType(ss);
            properties.bondIndex = util::extractFromStream(ss, int());
            properties.angleIndex = util::extractFromStream(ss, int());
            properties.dihedralIndex = util::extractFromStream(ss, int());
            properties.bondLength = util::extractFromStream(ss, double());
            properties.angle = util::extractFromStream(ss, double());
            properties.dihedral = util::extractFromStream(ss, double());
            atom->setCharge(util::extractFromStream(ss, double()));
        }

        void determine3dCoordinate(Atom* atom, const PrepAtomProperties& properties)
        {
            std::vector<Atom*> foundAtoms = findDihedralAtoms(atom);
            atom->setCoordinate(calculateCoordinateFromInternalCoords(
                foundAtoms.at(3)->coordinate(),
                foundAtoms.at(2)->coordinate(),
                foundAtoms.at(1)->coordinate(),
                properties.angle,
                properties.dihedral,
                properties.bondLength));
        }

        void print(Atom* atom, const PrepAtomProperties& properties, std::ostream& out)
        {
            out << std::setw(3) << atom->getNumber() << std::setw(6) << atom->getName() << std::setw(6)
                << atom->getType() << std::setw(3) << topologicalTypeNames[properties.topologicalType] << std::setw(4)
                << properties.bondIndex << std::setw(4) << properties.angleIndex << std::setw(4)
                << properties.dihedralIndex << std::setw(10) << properties.bondLength << std::setw(10)
                << properties.angle << std::setw(10) << properties.dihedral << std::setw(10) << atom->getCharge();
        }

        void write(Atom* atom, const PrepAtomProperties& properties, std::ostream& stream)
        {
            stream << std::right << std::setw(2) << atom->getNumber() << " " << std::left << std::setw(4)
                   << atom->getName() << " " << std::left << std::setw(3) << atom->getType() << " " << std::setw(1)
                   << topologicalTypeNames[properties.topologicalType] << " " << std::right << std::setw(2)
                   << properties.bondIndex << " " << std::right << std::setw(2) << properties.angleIndex << " "
                   << std::right << std::setw(2) << properties.dihedralIndex << " ";
            stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.bondLength << " ";
            stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.angle << " ";
            stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.dihedral;
            stream << "    " << std::right << std::setw(8) << std::fixed << std::setprecision(4) << atom->getCharge()
                   << std::endl;
        }
    } // namespace prep
} // namespace gmml
