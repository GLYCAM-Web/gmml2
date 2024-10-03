#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINOVERLAPRESOLUTION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINOVERLAPRESOLUTION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>

struct GlycoproteinState
{
    cds::Overlap overlap;
    std::vector<size_t> overlapSites;
    std::vector<std::vector<cds::ResidueLinkageShapePreference>> preferences;
    std::vector<std::vector<std::vector<cds::AngleWithMetadata>>> shape;
};

struct OverlapWeight
{
    double protein;
    double glycan;
    double self;
};

struct OverlapResidues
{
    std::vector<cds::Residue*> protein;
    std::vector<cds::Residue*> glycan;
};

typedef std::function<std::vector<std::vector<cds::AngleWithMetadata>>(
    const AssemblyGraphs&, AssemblyData&, size_t glycanId, OverlapWeight,
    std::vector<cds::ResidueLinkageShapePreference>&)>
    WiggleGlycan;

typedef std::function<std::vector<cds::ResidueLinkageShapePreference>(const std::vector<cds::ResidueLinkage>&)>
    LinkageShapeRandomizer;

void setLinkageShape(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                     const std::vector<std::vector<std::vector<cds::AngleWithMetadata>>>& shape);
void setLinkageShapeToPreference(const AssemblyGraphs& graphs, AssemblyData& data, size_t linkageId,
                                 const cds::ResidueLinkageShapePreference& preference);
cds::Overlap intraGlycanOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t glycanId);
cds::Overlap moleculeOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t moleculeA,
                              size_t moleculeB);
cds::Overlap moleculeResidueOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t molecule,
                                     size_t residue);
cds::Overlap totalOverlaps(OverlapWeight weight, const AssemblyGraphs& graphs, const AssemblyData& data);
std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const AssemblyGraphs& graphs,
                                              const AssemblyData& data);
std::vector<cds::AngleWithMetadata> wiggleLinkage(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                                                  size_t linkageId, const cds::AngleSearchSettings& searchSettings,
                                                  OverlapWeight weight,
                                                  const cds::ResidueLinkageShapePreference& shapePreference);
std::vector<std::vector<cds::AngleWithMetadata>>
wiggleGlycan(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
             const cds::AngleSearchSettings& searchSettings, OverlapWeight weight,
             const std::vector<cds::ResidueLinkageShapePreference>& preferences);
void updateGlycanBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId);
GlycoproteinState randomDescent(pcg32 rng, LinkageShapeRandomizer randomizeShape, WiggleGlycan wiggleGlycan,
                                int persistCycles, OverlapWeight overlapWeight, const AssemblyGraphs& graphs,
                                AssemblyData& data, std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                const GlycoproteinState& initialState);
#endif
