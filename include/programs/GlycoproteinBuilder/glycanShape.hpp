#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GLYCANSHAPE_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GLYCANSHAPE_HPP

#include "include/CentralDataStructure/Shapers/dihedralAngleSearchTypes.hpp"
#include "include/CentralDataStructure/Shapers/dihedralShapeTypes.hpp"
#include "include/External_Libraries/PCG/pcg_random.h"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/programs/GlycoproteinBuilder/glycoproteinStructs.hpp"

namespace gmml
{
    namespace gpbuilder
    {
        std::array<Coordinate, 4> dihedralCoordinates(
            const AssemblyData& data, const assembly::Bounds& bounds, size_t dihedralId);
        AngleWithMetadata currentDihedralShape(
            const AssemblyData& data, const MutableData& mutableData, size_t dihedralId);

        void setDihedralAngle(
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t dihedralId,
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
