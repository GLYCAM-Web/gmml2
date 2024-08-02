#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPE_HPP

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <functional>

namespace cds
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

    typedef std::function<size_t(DihedralAngleDataVector)> MetadataDistribution;
    typedef std::function<double(DihedralAngleData)> AngleDistribution;

    void setDihedralAngle(RotatableDihedral& dihedral, cds::AngleWithMetadata target);
    bool setSpecificShape(RotatableDihedral& dihedral, const DihedralAngleDataVector& metadataVector,
                          std::string dihedralName, std::string selectedRotamer);
    std::vector<AngleWithMetadata> currentShape(const std::vector<RotatableDihedral>& dihedrals,
                                                const std::vector<DihedralAngleDataVector>& metadata);
    std::vector<std::vector<AngleWithMetadata>> currentShape(const std::vector<ResidueLinkage>& linkages);
    void setShape(std::vector<ResidueLinkage>& linkages, const std::vector<std::vector<AngleWithMetadata>>& angles);

    void setShape(std::vector<RotatableDihedral>& dihedrals, const std::vector<AngleWithMetadata>& angles);
    void setDefaultShapeUsingMetadata(std::vector<cds::RotatableDihedral>& dihedrals,
                                      const std::vector<DihedralAngleDataVector>& metadata);
    void setRandomShapeUsingMetadata(MetadataDistribution randomMetadata, AngleDistribution randomAngle,
                                     gmml::MolecularMetadata::GLYCAM::RotamerType rotamerType,
                                     std::vector<RotatableDihedral>& dihedrals,
                                     const std::vector<DihedralAngleDataVector>& metadata);
    void setRandomShapeUsingMetadata(MetadataDistribution randomMetadata, AngleDistribution randomAngle,
                                     std::vector<ResidueLinkage>& linkages);
    void setSpecificShapeUsingMetadata(std::vector<RotatableDihedral>& dihedrals,
                                       const std::vector<DihedralAngleDataVector>& metadata, size_t shapeNumber);
    void setSpecificShape(std::vector<RotatableDihedral>& dihedrals,
                          const std::vector<DihedralAngleDataVector>& metadata, std::string dihedralName,
                          std::string selectedRotamer);

} // namespace cds
#endif
