#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <sstream>

using cds::RotatableDihedral;

cds::DihedralCoordinates cds::dihedralCoordinates(const cds::RotatableDihedral& dihedral)
{
    auto& atoms                       = dihedral.atoms;
    std::array<Coordinate*, 4> coords = {atoms[3]->getCoordinate(), atoms[2]->getCoordinate(),
                                         atoms[1]->getCoordinate(), atoms[0]->getCoordinate()};
    return DihedralCoordinates {*coords[0], *coords[1], *coords[2], *coords[3]};
}

void cds::setDihedralAngle(RotatableDihedral& dihedral, cds::AngleWithMetadata target)
{
    auto matrix = rotationTo(dihedralCoordinates(dihedral), constants::degree2Radian(target.value));
    matrix.rotateCoordinates(dihedral.movingCoordinates);
    dihedral.currentMetadata = target.metadata;
}

bool cds::setSpecificShape(RotatableDihedral& dihedral, std::string dihedralName, std::string selectedRotamer)
{
    auto& metadataVector = dihedral.metadataVector;
    if (dihedralName == metadataVector[0].dihedral_angle_name_)
    {
        for (auto& metadata : metadataVector)
        {
            if (metadata.rotamer_name_ == selectedRotamer)
            {
                setDihedralAngle(dihedral, {metadata.default_angle_value_, metadata});
                return true;
            }
        }
    }
    return false;
}

std::string cds::print(const RotatableDihedral& dihedral)
{
    auto& atoms = dihedral.atoms;
    std::stringstream ss;
    ss << atoms[0]->getName() << ", " << atoms[1]->getName() << ", " << atoms[2]->getName() << ", "
       << atoms[3]->getName() << ": " << cds::angle(dihedralCoordinates(dihedral)) << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
}
