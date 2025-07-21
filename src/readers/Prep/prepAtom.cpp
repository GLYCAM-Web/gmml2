#include "include/readers/Prep/prepAtom.hpp"

#include "include/geometry/measurements.hpp" //get_cartesian_point_from_internal_coords()
#include "include/util/casting.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <iomanip>
#include <ios>
#include <sstream>

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
                if (s == "M")
                {
                    return kTopTypeM;
                }
                else if (s == "S")
                {
                    return kTopTypeS;
                }
                else if (s == "B")
                {
                    return kTopTypeB;
                }
                else if (s == "E")
                {
                    return kTopTypeE;
                }
                else
                {
                    return kTopType3;
                }
            }

            std::string getStringFormatOfTopologicalType(TopologicalType type)
            {
                switch (type)
                {
                    case kTopTypeE:
                        return "E";
                    case kTopTypeS:
                        return "S";
                    case kTopTypeB:
                        return "B";
                    case kTopType3:
                        return "3";
                    case kTopType4:
                        return "4";
                    case kTopTypeM:
                        return "M";
                    default:
                        return "E";
                }
            }

            std::vector<PrepAtom*> findDihedralAtoms(PrepAtom* initial)
            {
                std::vector<PrepAtom*> result = {initial};
                for (int currentDepth = 0; currentDepth < 3; currentDepth++)
                {
                    PrepAtom* parent = util::erratic_cast<PrepAtom*>(
                        result.back()->getParents().front()); // Go up the first parent only. Loops may create another
                                                              // parent, but they should be ignored.
                    result.push_back(parent);
                }
                return result;
            }
        } // namespace

        //////////////////////////////////////////////////////////
        //                       Constructor                    //
        //////////////////////////////////////////////////////////
        PrepAtom::PrepAtom(const std::string& line)
        {
            std::stringstream ss(line);
            this->setNumber(util::extractFromStream(ss, int()));
            this->setName(util::extractFromStream(ss, std::string()));
            this->setType(util::extractFromStream(ss, std::string()));
            properties.topologicalType = extractAtomTopologicalType(ss);
            properties.bondIndex = util::extractFromStream(ss, int());
            properties.angleIndex = util::extractFromStream(ss, int());
            properties.dihedralIndex = util::extractFromStream(ss, int());
            properties.bondLength = util::extractFromStream(ss, double());
            properties.angle = util::extractFromStream(ss, double());
            properties.dihedral = util::extractFromStream(ss, double());
            this->setCharge(util::extractFromStream(ss, double()));
        }

        // Move Ctor
        PrepAtom::PrepAtom(PrepAtom&& other) noexcept : PrepAtom() //: Atom(other)
        {
            swap(*this, other);
        }

        PrepAtom::PrepAtom(const PrepAtom& other) noexcept : Atom(other), properties(other.properties) {}

        PrepAtom& PrepAtom::operator=(PrepAtom other) noexcept
        {
            swap(*this, other);
            return *this;
        }

        void PrepAtom::Determine3dCoordinate()
        {
            // std::cout << "Determining 3d Coordinates for " << this->getName() << "\n";
            std::vector<PrepAtom*> foundAtoms = findDihedralAtoms(this);
            this->setCoordinate(calculateCoordinateFromInternalCoords(
                foundAtoms.at(3)->coordinate(),
                foundAtoms.at(2)->coordinate(),
                foundAtoms.at(1)->coordinate(),
                properties.angle,
                properties.dihedral,
                properties.bondLength));
            return;
        }

        //////////////////////////////////////////////////////////
        //                     DISPLAY FUNCTIONS                //
        //////////////////////////////////////////////////////////
        void PrepAtom::Print(std::ostream& out) const
        {
            out << std::setw(3) << this->getNumber() << std::setw(6) << this->getName() << std::setw(6)
                << this->getType();
            if (properties.topologicalType == kTopTypeE)
            {
                out << std::setw(3) << "E";
            }
            else if (properties.topologicalType == kTopTypeS)
            {
                out << std::setw(3) << "S";
            }
            else if (properties.topologicalType == kTopTypeB)
            {
                out << std::setw(3) << "B";
            }
            else if (properties.topologicalType == kTopType3)
            {
                out << std::setw(3) << "3";
            }
            else if (properties.topologicalType == kTopType4)
            {
                out << std::setw(3) << "4";
            }
            else if (properties.topologicalType == kTopTypeM)
            {
                out << std::setw(3) << "M";
            }
            else
            {
                out << std::setw(3) << "-";
            }

            out << std::setw(4) << properties.bondIndex << std::setw(4) << properties.angleIndex << std::setw(4)
                << properties.dihedralIndex << std::setw(10) << properties.bondLength << std::setw(10)
                << properties.angle << std::setw(10) << properties.dihedral << std::setw(10) << this->getCharge();
            //        << endl;
        }

        void PrepAtom::Write(std::ostream& stream) const
        {
            stream << std::right << std::setw(2) << this->getNumber() << " " << std::left << std::setw(4)
                   << this->getName() << " " << std::left << std::setw(3) << this->getType() << " " << std::setw(1)
                   << getStringFormatOfTopologicalType(properties.topologicalType) << " " << std::right << std::setw(2)
                   << properties.bondIndex << " " << std::right << std::setw(2) << properties.angleIndex << " "
                   << std::right << std::setw(2) << properties.dihedralIndex << " ";
            stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.bondLength << " ";
            stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.angle << " ";
            stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.dihedral;
            stream << "    " << std::right << std::setw(8) << std::fixed << std::setprecision(4) << this->getCharge()
                   << std::endl;
        }
    } // namespace prep
} // namespace gmml
