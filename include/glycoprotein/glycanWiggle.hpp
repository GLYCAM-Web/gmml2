#ifndef INCLUDE_GLYCOPROTEIN_GLYCANWIGGLE_HPP
#define INCLUDE_GLYCOPROTEIN_GLYCANWIGGLE_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralAngleSearchTypes.hpp"
#include "include/carbohydrate/dihedralShapeTypes.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"

namespace gmml
{
    namespace gpbuilder
    {
        void wiggleLinkage(
            const DihedralAngleDataTable& dihedralAngleTable,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            MutableData& mutableData,
            size_t linkageId,
            const AngleSearchSettings& searchSettings,
            const ResidueLinkageShapePreference& shapePreference);

        void wiggleGlycan(
            const DihedralAngleDataTable& dihedralAngleTable,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const AngleSearchSettings& searchSettings,
            const GlycanShapePreference& preferences,
            MutableData& mutableData,
            size_t glycanId);
    } // namespace gpbuilder
} // namespace gmml

#endif
