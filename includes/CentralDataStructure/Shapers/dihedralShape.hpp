#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPE_HPP

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

// This class stores the four atoms that define a dihedral angle, the atoms that move when it is rotated
// and, if moved, the previous dihedral angle, which allows me to reset easily.
namespace cds
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

    void setDihedralAngle(RotatableDihedral& dihedral, AngleWithMetadata target);
    bool setSpecificShape(RotatableDihedral& dihedral, std::string dihedralName, std::string selectedRotamer);
    std::vector<AngleWithMetadata> currentShape(const std::vector<RotatableDihedral>& dihedrals);
    std::vector<std::vector<AngleWithMetadata>> currentShape(const std::vector<ResidueLinkage>& linkages);
    void setShape(std::vector<ResidueLinkage>& linkages, const std::vector<std::vector<AngleWithMetadata>>& angles);

    void setShape(std::vector<RotatableDihedral>& dihedrals, const std::vector<AngleWithMetadata>& angles);
    void setDefaultShapeUsingMetadata(std::vector<cds::RotatableDihedral>& dihedrals);
    void setRandomShapeUsingMetadata(std::vector<RotatableDihedral>& dihedrals);
    void setSpecificShapeUsingMetadata(std::vector<RotatableDihedral>& dihedrals, int shapeNumber);
    void setSpecificShape(std::vector<RotatableDihedral>& dihedrals, std::string dihedralName,
                          std::string selectedRotamer);

} // namespace cds
#endif
