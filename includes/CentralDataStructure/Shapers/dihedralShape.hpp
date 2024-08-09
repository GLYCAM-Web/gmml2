#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPE_HPP

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <vector>
#include <string>
#include <functional>
#include <variant>

namespace cds
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

    typedef std::function<std::vector<size_t>(DihedralAngleDataVector)> MetadataDistribution;
    typedef std::function<double(DihedralAngleData)> AngleDistribution;

    typedef std::function<std::vector<size_t>(std::vector<size_t>, const RotatableDihedral&)>
        MetadataPreferenceSelection;

    struct ConformerShapePreference
    {
        std::vector<std::vector<double>> angles;
        std::vector<size_t> metadataOrder;
    };

    struct PermutationShapePreference
    {
        std::vector<std::vector<double>> angles;
        std::vector<std::vector<size_t>> metadataOrder;
    };

    typedef std::variant<ConformerShapePreference, PermutationShapePreference> ResidueLinkageShapePreference;

    void setDihedralAngle(RotatableDihedral& dihedral, cds::AngleWithMetadata target);
    bool setSpecificShape(RotatableDihedral& dihedral, const DihedralAngleDataVector& metadataVector,
                          std::string dihedralName, std::string selectedRotamer);
    void setSpecificShape(std::vector<RotatableDihedral>& dihedrals,
                          const std::vector<DihedralAngleDataVector>& metadata, std::string dihedralName,
                          std::string selectedRotamer);
    std::vector<AngleWithMetadata> currentShape(const std::vector<RotatableDihedral>& dihedrals,
                                                const std::vector<DihedralAngleDataVector>& metadata);
    std::vector<std::vector<AngleWithMetadata>> currentShape(const std::vector<ResidueLinkage>& linkages);
    void setShape(std::vector<ResidueLinkage>& linkages, const std::vector<std::vector<AngleWithMetadata>>& angles);
    void setShape(std::vector<RotatableDihedral>& dihedrals, const std::vector<AngleWithMetadata>& angles);
    void setShapeToPreference(ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference);
    void setShapeToPreference(std::vector<ResidueLinkage>& linkages,
                              const std::vector<ResidueLinkageShapePreference>& preferences);
    ResidueLinkageShapePreference linkageShapePreference(MetadataDistribution metadataDistribution,
                                                         AngleDistribution angleDistribution,
                                                         const ResidueLinkage& linkage);
    ResidueLinkageShapePreference defaultShapePreference(const ResidueLinkage& linkage);
    ResidueLinkageShapePreference selectedRotamersOnly(MetadataPreferenceSelection metadataSelection,
                                                       const ResidueLinkage& linkage,
                                                       const ResidueLinkageShapePreference& preference);
    ResidueLinkageShapePreference firstRotamerOnly(const ResidueLinkage& linkage,
                                                   const ResidueLinkageShapePreference& preference);
    ResidueLinkageShapePreference currentRotamerOnly(const ResidueLinkage& linkage,
                                                     const ResidueLinkageShapePreference& preference);
    std::vector<ResidueLinkageShapePreference> linkageShapePreference(MetadataDistribution metadataDistribution,
                                                                      AngleDistribution angleDistribution,
                                                                      const std::vector<ResidueLinkage>& linkages);
    std::vector<ResidueLinkageShapePreference> defaultShapePreference(const std::vector<ResidueLinkage>& linkages);
    std::vector<ResidueLinkageShapePreference>
    firstRotamerOnly(const std::vector<ResidueLinkage>& linkages,
                     const std::vector<ResidueLinkageShapePreference>& preferences);
    std::vector<ResidueLinkageShapePreference>
    currentRotamerOnly(const std::vector<ResidueLinkage>& linkages,
                       const std::vector<ResidueLinkageShapePreference>& preferences);

} // namespace cds
#endif
