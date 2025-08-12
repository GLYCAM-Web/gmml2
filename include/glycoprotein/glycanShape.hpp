#ifndef INCLUDE_GLYCOPROTEIN_GLYCANSHAPE_HPP
#define INCLUDE_GLYCOPROTEIN_GLYCANSHAPE_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralAngleSearchTypes.hpp"
#include "include/carbohydrate/dihedralShapeTypes.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/geometry/overlap.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"

namespace gmml
{
    namespace gpbuilder
    {
        std::array<Coordinate, 4> dihedralCoordinates(
            const AssemblyData& data, const assembly::Bounds& bounds, size_t bondId);
        AngleWithMetadata currentDihedralShape(
            const DihedralAngleDataTable& dihedralAngleTable,
            const AssemblyData& data,
            const MutableData& mutableData,
            size_t bondId);

        void setDihedralAngle(
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t bondId,
            const AngleWithMetadata& target);

        void setLinkageShape(
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t glycanId,
            const std::vector<AngleWithMetadata>& recordedShape);

        void setLinkageShapeToPreference(
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t linkageId,
            const ResidueLinkageShapePreference& preference);
    } // namespace gpbuilder
} // namespace gmml

#endif
