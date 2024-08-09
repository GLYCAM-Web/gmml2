#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINOVERLAPRESOLUTION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINOVERLAPRESOLUTION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>

struct ResiduesByType
{
    std::vector<std::vector<cds::Residue*>> protein;
    std::vector<std::vector<cds::Residue*>> glycan;
    std::vector<std::vector<cds::Residue*>> all;
};

typedef std::function<std::vector<cds::ResidueLinkageShapePreference>(std::vector<cds::ResidueLinkage>)>
    LinkageShapeRandomizer;

cds::Overlap countGlycositeOverlaps(const std::vector<Residue*>& overlapResidues,
                                    const std::vector<Residue*>& glycositeResidues);
cds::Overlap countTotalOverlaps(const std::vector<std::vector<Residue*>>& overlapResidues,
                                const std::vector<std::vector<Residue*>>& glycositeResidues);
std::vector<size_t> determineSitesWithOverlap(uint overlapTolerance,
                                              const std::vector<std::vector<Residue*>>& overlapResidues,
                                              const std::vector<std::vector<Residue*>>& glycositeResidues);
cds::Overlap
wiggleSitesWithOverlaps(pcg32& rng, uint overlapTolerance, int persistCycles, bool firstLinkageOnly, int interval,
                        std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                        const std::vector<std::vector<cds::ResidueLinkageShapePreference>>& glycositePreferences,
                        const std::vector<std::vector<Residue*>>& overlapResidues,
                        const std::vector<std::vector<Residue*>>& glycositeResidues);
cds::Overlap randomDescent(pcg32& rng, LinkageShapeRandomizer randomizeShape, bool monte_carlo, int persistCycles,
                           uint overlapTolerance, std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                           std::vector<std::vector<cds::ResidueLinkageShapePreference>>& glycositePreferences,
                           const ResiduesByType& overlapResidues,
                           const std::vector<std::vector<Residue*>>& glycositeResidues);
bool dumbRandomWalk(LinkageShapeRandomizer randomizeShape, uint overlapTolerance, int maxCycles,
                    std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                    std::vector<std::vector<cds::ResidueLinkageShapePreference>>& glycositePreferences,
                    const std::vector<std::vector<Residue*>>& overlapResidues,
                    const std::vector<std::vector<Residue*>>& glycositeResidues);
#endif
