#ifndef INCLUDE_CARBOHYDRATE_DIHEDRALSHAPE_HPP
#define INCLUDE_CARBOHYDRATE_DIHEDRALSHAPE_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageTypes.hpp"
#include "include/carbohydrate/dihedralShapeTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/orientation.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <functional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    template<typename T>
    T onResidueLinkageShapePreference(
        std::function<T(const ConformerShapePreference&)>& onConformer,
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
            throw std::runtime_error("unhandled linkage shape preference in onResidueLinkageShapePreference");
        }
    };

    void setDihedralAngle(RotatableBond& bond, AngleWithMetadata target);

    bool setSpecificShape(
        const DihedralAngleDataTable& metadataTable,
        RotatableBond& dihedral,
        const std::vector<size_t>& metadataVector,
        std::string dihedralName,
        std::string selectedRotamer);

    void setSpecificShape(
        const DihedralAngleDataTable& metadataTable,
        std::vector<RotatableBond>& bonds,
        const std::vector<std::vector<size_t>>& metadata,
        std::string dihedralName,
        std::string selectedRotamer);

    std::vector<AngleWithMetadata> currentShape(
        const DihedralAngleDataTable& metadataTable,
        const std::vector<RotatableBond>& bonds,
        const std::vector<std::vector<size_t>>& metadata);

    std::vector<std::vector<AngleWithMetadata>> currentShape(
        const DihedralAngleDataTable& metadataTable, const std::vector<ResidueLinkage>& linkages);
    void setShape(std::vector<ResidueLinkage>& linkages, const std::vector<std::vector<AngleWithMetadata>>& angles);
    void setShape(std::vector<RotatableBond>& bonds, const std::vector<AngleWithMetadata>& angles);
    void setShapeToPreference(ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference);
    void setShapeToPreference(std::vector<ResidueLinkage>& linkages, const GlycanShapePreference& preferences);

    ResidueLinkageShapePreference linkageShapePreference(
        MetadataDistribution metadataDistribution,
        AngleDistribution angleDistribution,
        const DihedralAngleDataTable& metadataTable,
        const RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& dihedralMetadata);

    ResidueLinkageShapePreference defaultShapePreference(
        const DihedralAngleDataTable& metadataTable,
        const RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& dihedralMetadata);

    ResidueLinkageShapePreference selectedRotamersOnly(
        MetadataPreferenceSelection metadataSelection,
        const ResidueLinkage& linkage,
        const ResidueLinkageShapePreference& preference);

    ResidueLinkageShapePreference firstRotamerOnly(
        const ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference);

    ResidueLinkageShapePreference currentRotamerOnly(
        const ResidueLinkage& linkage, const ResidueLinkageShapePreference& preference);

    GlycanShapePreference linkageShapePreference(
        MetadataDistribution metadataDistribution,
        AngleDistribution angleDistribution,
        const DihedralAngleDataTable& metadataTable,
        const std::vector<ResidueLinkage>& linkages);

    GlycanShapePreference defaultShapePreference(
        const DihedralAngleDataTable& metadataTable, const std::vector<ResidueLinkage>& linkages);
    GlycanShapePreference firstRotamerOnly(
        const std::vector<ResidueLinkage>& linkages, const GlycanShapePreference& preferences);
    GlycanShapePreference currentRotamerOnly(
        const std::vector<ResidueLinkage>& linkages, const GlycanShapePreference& preferences);

} // namespace gmml
#endif
