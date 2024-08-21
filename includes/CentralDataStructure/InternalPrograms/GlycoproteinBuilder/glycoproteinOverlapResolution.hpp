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

typedef std::function<std::vector<cds::ResidueLinkageShapePreference>(std::vector<cds::ResidueLinkage>)>
    LinkageShapeRandomizer;

cds::Overlap countTotalOverlaps(const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues);
std::vector<size_t> determineSitesWithOverlap(uint overlapTolerance,
                                              const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                              const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues);
GlycoproteinState wiggleSitesWithOverlaps(pcg32& rng, cds::SearchAngles searchAngles, uint overlapTolerance,
                                          int persistCycles, bool firstLinkageOnly,
                                          std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                          const GlycoproteinState& initialState,
                                          const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                          const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues);
GlycoproteinState randomDescent(pcg32& rng, LinkageShapeRandomizer randomizeShape, cds::SearchAngles searchAngles,
                                bool monte_carlo, int persistCycles, uint overlapTolerance,
                                std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                const GlycoproteinState& initialState,
                                const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues);
#endif
