#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCANSHAPE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCANSHAPE_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

namespace glycoproteinBuilder
{
    std::array<cds::Coordinate, 4> dihedralCoordinates(const AssemblyData& data, const assembly::Bounds& bounds,
                                                       size_t dihedralId);
    void setDihedralAngle(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                          size_t dihedralId, const cds::AngleWithMetadata& target);
    void setLinkageShape(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                         size_t glycanId, const std::vector<cds::AngleWithMetadata>& recordedShape);
    void setLinkageShapeToPreference(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                                     size_t linkageId, const cds::ResidueLinkageShapePreference& preference);
} // namespace glycoproteinBuilder
#endif
