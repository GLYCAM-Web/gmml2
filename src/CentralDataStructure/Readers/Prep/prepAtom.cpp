#include "includes/CentralDataStructure/Readers/Prep/prepAtom.hpp"

#include "includes/CentralDataStructure/Measurements/measurements.hpp" //get_cartesian_point_from_internal_coords()
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <iomanip>
#include <ios>
#include <sstream>

namespace
{
    prep::TopologicalType extractAtomTopologicalType(std::istream& ss)
    {
        std::string s;
        ss >> s;
        if (s == "M")
        {
            return prep::kTopTypeM;
        }
        else if (s == "S")
        {
            return prep::kTopTypeS;
        }
        else if (s == "B")
        {
            return prep::kTopTypeB;
        }
        else if (s == "E")
        {
            return prep::kTopTypeE;
        }
        else
        {
            return prep::kTopType3;
        }
    }

    std::string getStringFormatOfTopologicalType(prep::TopologicalType type)
    {
        switch (type)
        {
            case prep::kTopTypeE:
                return "E";
            case prep::kTopTypeS:
                return "S";
            case prep::kTopTypeB:
                return "B";
            case prep::kTopType3:
                return "3";
            case prep::kTopType4:
                return "4";
            case prep::kTopTypeM:
                return "M";
            default:
                return "E";
        }
    }

    std::vector<prep::PrepAtom*> findDihedralAtoms(prep::PrepAtom* initial)
    {
        std::vector<prep::PrepAtom*> result = {initial};
        for (int currentDepth = 0; currentDepth < 3; currentDepth++)
        {
            prep::PrepAtom* parent = codeUtils::erratic_cast<prep::PrepAtom*>(
                result.back()->getParents().front()); // Go up the first parent only. Loops may create another parent,
                                                      // but they should be ignored.
            result.push_back(parent);
        }
        return result;
    }
} // namespace

using prep::PrepAtom;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepAtom::PrepAtom(const std::string& line)
{
    std::stringstream ss(line);
    this->setNumber(codeUtils::extractFromStream(ss, int()));
    this->setName(codeUtils::extractFromStream(ss, std::string()));
    this->setType(codeUtils::extractFromStream(ss, std::string()));
    properties.topologicalType = extractAtomTopologicalType(ss);
    properties.bondIndex = codeUtils::extractFromStream(ss, int());
    properties.angleIndex = codeUtils::extractFromStream(ss, int());
    properties.dihedralIndex = codeUtils::extractFromStream(ss, int());
    properties.bondLength = codeUtils::extractFromStream(ss, double());
    properties.angle = codeUtils::extractFromStream(ss, double());
    properties.dihedral = codeUtils::extractFromStream(ss, double());
    this->setCharge(codeUtils::extractFromStream(ss, double()));
}

// Move Ctor
PrepAtom::PrepAtom(PrepAtom&& other) noexcept : PrepAtom() //: cds::Atom(other)
{
    swap(*this, other);
}

PrepAtom::PrepAtom(const PrepAtom& other) noexcept : cds::Atom(other), properties(other.properties) {}

PrepAtom& PrepAtom::operator=(PrepAtom other) noexcept
{
    swap(*this, other);
    return *this;
}

void PrepAtom::Determine3dCoordinate()
{
    // std::cout << "Determining 3d Coordinates for " << this->getName() << "\n";
    std::vector<PrepAtom*> foundAtoms = findDihedralAtoms(this);
    this->setCoordinate(cds::calculateCoordinateFromInternalCoords(
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
    out << std::setw(3) << this->getNumber() << std::setw(6) << this->getName() << std::setw(6) << this->getType();
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
        << properties.dihedralIndex << std::setw(10) << properties.bondLength << std::setw(10) << properties.angle
        << std::setw(10) << properties.dihedral << std::setw(10) << this->getCharge();
    //        << endl;
}

void PrepAtom::Write(std::ostream& stream) const
{
    stream << std::right << std::setw(2) << this->getNumber() << " " << std::left << std::setw(4) << this->getName()
           << " " << std::left << std::setw(3) << this->getType() << " " << std::setw(1)
           << getStringFormatOfTopologicalType(properties.topologicalType) << " " << std::right << std::setw(2)
           << properties.bondIndex << " " << std::right << std::setw(2) << properties.angleIndex << " " << std::right
           << std::setw(2) << properties.dihedralIndex << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.bondLength << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.angle << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << properties.dihedral;
    stream << "    " << std::right << std::setw(8) << std::fixed << std::setprecision(4) << this->getCharge()
           << std::endl;
}
