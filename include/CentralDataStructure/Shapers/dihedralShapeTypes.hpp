#ifndef INCLUDE_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPETYPES_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALSHAPETYPES_HPP

#include "include/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <functional>
#include <string>
#include <variant>
#include <vector>

namespace gmml
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

    typedef std::function<std::vector<size_t>(const DihedralAngleDataTable&, const std::vector<size_t>)>
        MetadataDistribution;
    typedef std::function<double(const DihedralAngleData&)> AngleDistribution;

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

    typedef std::vector<ResidueLinkageShapePreference> GlycanShapePreference;
} // namespace gmml
#endif
