#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCANSHAPE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCANSHAPE_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

namespace glycoproteinBuilder
{
    void updateResidueBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index);
    void updateResidueMoleculeBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index);
    void updateMoleculeBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index);
    void updateGlycanBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId);
    std::array<cds::Coordinate, 4> dihedralCoordinates(const AssemblyGraphs& graphs, const AssemblyData& data,
                                                       size_t dihedralId);
    void setDihedralAngle(const AssemblyGraphs& graphs, AssemblyData& data, size_t linkageId, size_t dihedralId,
                          const cds::AngleWithMetadata& target);
    void setLinkageShape(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                         const std::vector<cds::AngleWithMetadata>& recordedShape);
    void setLinkageShapeToPreference(const AssemblyGraphs& graphs, AssemblyData& data, size_t linkageId,
                                     const cds::ResidueLinkageShapePreference& preference);
} // namespace glycoproteinBuilder
#endif
