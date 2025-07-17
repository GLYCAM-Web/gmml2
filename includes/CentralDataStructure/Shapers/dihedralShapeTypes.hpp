#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPETYPES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPETYPES_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
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
} // namespace cds
#endif
