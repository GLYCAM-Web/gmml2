#ifndef INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCANWIGGLE_HPP
#define INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCANWIGGLE_HPP

#include "include/CentralDataStructure/Shapers/dihedralAngleSearchTypes.hpp"
#include "include/CentralDataStructure/Shapers/dihedralShapeTypes.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"

namespace gmml
{
    namespace gpbuilder
    {
        void wiggleLinkage(
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            MutableData& mutableData,
            size_t linkageId,
            const AngleSearchSettings& searchSettings,
            const ResidueLinkageShapePreference& shapePreference);

        void wiggleGlycan(
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
