#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINOVERLAPRESOLUTION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINOVERLAPRESOLUTION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>

struct GlycoproteinState
{
    cds::Overlap overlap;
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

typedef std::function<std::vector<cds::ResidueLinkageShapePreference>(std::vector<cds::ResidueLinkage>)>
    LinkageShapeRandomizer;
cds::Overlap intraGlycanOverlaps(const std::vector<cds::ResidueLinkage>& linkages);
cds::Overlap countOverlaps(const std::vector<Residue*>& overlapResidues,
                           const cds::ResiduesWithOverlapWeight& glycositeResidues);
cds::Overlap siteOverlaps(OverlapWeight weight, const OverlapResidues& overlapResidues,
                          const cds::ResiduesWithOverlapWeight& glycositeResidues,
                          const std::vector<cds::ResidueLinkage>& linkages);
cds::Overlap totalOverlaps(OverlapWeight weight, const std::vector<OverlapResidues>& overlapResidues,
                           const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues,
                           const std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages);
std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites,
                                              const std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                              const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                              const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues);
std::vector<cds::AngleWithMetadata> wiggleLinkage(cds::SearchAngles searchAngles, cds::ResidueLinkage& linkage,
                                                  const cds::ResidueLinkageShapePreference& shapePreference,
                                                  const std::array<cds::ResiduesWithOverlapWeight, 2> overlapInput);
std::vector<std::vector<cds::AngleWithMetadata>>
wiggleGlycosite(cds::SearchAngles searchAngles, OverlapWeight weight, std::vector<cds::ResidueLinkage>& linkages,
                const std::vector<cds::ResidueLinkageShapePreference>& preferences,
                const OverlapResidues& overlapResidues);
GlycoproteinState randomDescent(pcg32 rng, LinkageShapeRandomizer randomizeShape, cds::SearchAngles searchAngles,
                                std::function<bool(int)> acceptViaMetropolisCriterion, int persistCycles,
                                OverlapWeight overlapWeight,
                                std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                const GlycoproteinState& initialState,
                                const std::vector<OverlapResidues>& overlapResidues,
                                const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues);
#endif
