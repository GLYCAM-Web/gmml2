#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPE_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <vector>
#include <string>
#include <functional>
#include <variant>

namespace cds
{
    struct Bounds
    {
        double lower;
        double upper;
    };

    struct AngleWithMetadata
    {
        double value;
        double preference;
        size_t metadataIndex;
    };

    typedef std::function<std::vector<size_t>(const GlycamMetadata::DihedralAngleDataTable&, const std::vector<size_t>)>
        MetadataDistribution;
    typedef std::function<double(const GlycamMetadata::DihedralAngleData&)> AngleDistribution;

    typedef std::function<std::vector<size_t>(std::vector<size_t>, const RotatableDihedral&)>
        MetadataPreferenceSelection;

    struct ConformerShapePreference
    {
        std::vector<bool> isFrozen;
        std::vector<std::vector<double>> angles;
        std::vector<size_t> metadataOrder;
    };

    struct PermutationShapePreference
    {
        std::vector<std::vector<double>> angles;
        std::vector<std::vector<size_t>> metadataOrder;
    };

    typedef std::variant<ConformerShapePreference, PermutationShapePreference> ResidueLinkageShapePreference;

    typedef std::vector<cds::ResidueLinkageShapePreference> GlycanShapePreference;

    template<typename T>
    T onResidueLinkageShapePreference(std::function<T(const ConformerShapePreference&)>& onConformer,
                                      std::function<T(const PermutationShapePreference&)>& onPermutation,
                                      const ResidueLinkageShapePreference& preference)
    {
        if (std::holds_alternative<ConformerShapePreference>(preference))
        {
            return onConformer(std::get<ConformerShapePreference>(preference));
        }
        else if (std::holds_alternative<PermutationShapePreference>(preference))
        {
            return onPermutation(std::get<PermutationShapePreference>(preference));
        }
        else
        {
            throw std::runtime_error("unhandled linkage shape preference in cds::onResidueLinkageShapePreference");
        }
    };

    void setDihedralAngle(RotatableDihedral& dihedral, cds::AngleWithMetadata target);
    bool setSpecificShape(const GlycamMetadata::DihedralAngleDataTable& metadataTable, RotatableDihedral& dihedral,
                          const std::vector<size_t>& metadataVector, std::string dihedralName,
                          std::string selectedRotamer);
    void setSpecificShape(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                          std::vector<RotatableDihedral>& dihedrals, const std::vector<std::vector<size_t>>& metadata,
                          std::string dihedralName, std::string selectedRotamer);
    std::vector<AngleWithMetadata> currentShape(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                                const std::vector<RotatableDihedral>& dihedrals,
                                                const std::vector<std::vector<size_t>>& metadata);
    std::vector<std::vector<AngleWithMetadata>>
    currentShape(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                 const std::vector<ResidueLinkage>& linkages);
    void setShape(std::vector<ResidueLinkage>& linkages, const std::vector<std::vector<AngleWithMetadata>>& angles);
    void setShape(std::vector<RotatableDihedral>& dihedrals, const std::vector<AngleWithMetadata>& angles);
    void setShapeToPreference(ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference);
    void setShapeToPreference(std::vector<ResidueLinkage>& linkages, const GlycanShapePreference& preferences);
    ResidueLinkageShapePreference linkageShapePreference(MetadataDistribution metadataDistribution,
                                                         AngleDistribution angleDistribution,
                                                         const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                                         const GlycamMetadata::RotamerType rotamerType,
                                                         const std::vector<std::vector<size_t>>& dihedralMetadata);
    ResidueLinkageShapePreference defaultShapePreference(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                                         const GlycamMetadata::RotamerType rotamerType,
                                                         const std::vector<std::vector<size_t>>& dihedralMetadata);
    ResidueLinkageShapePreference selectedRotamersOnly(MetadataPreferenceSelection metadataSelection,
                                                       const ResidueLinkage& linkage,
                                                       const ResidueLinkageShapePreference& preference);
    ResidueLinkageShapePreference firstRotamerOnly(const ResidueLinkage& linkage,
                                                   const ResidueLinkageShapePreference& preference);
    ResidueLinkageShapePreference currentRotamerOnly(const ResidueLinkage& linkage,
                                                     const ResidueLinkageShapePreference& preference);
    GlycanShapePreference linkageShapePreference(MetadataDistribution metadataDistribution,
                                                 AngleDistribution angleDistribution,
                                                 const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                                 const std::vector<ResidueLinkage>& linkages);
    GlycanShapePreference defaultShapePreference(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                                 const std::vector<ResidueLinkage>& linkages);
    GlycanShapePreference firstRotamerOnly(const std::vector<ResidueLinkage>& linkages,
                                           const GlycanShapePreference& preferences);
    GlycanShapePreference currentRotamerOnly(const std::vector<ResidueLinkage>& linkages,
                                             const GlycanShapePreference& preferences);

} // namespace cds
#endif
