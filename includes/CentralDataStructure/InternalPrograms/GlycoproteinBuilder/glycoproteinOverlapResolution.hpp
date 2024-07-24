#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINOVERLAPRESOLUTION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINOVERLAPRESOLUTION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CentralDataStructure/residue.hpp"

struct ResiduesByType
{
    std::vector<std::vector<cds::Residue*>> protein;
    std::vector<std::vector<cds::Residue*>> glycan;
    std::vector<std::vector<cds::Residue*>> all;
};

cds::Overlap countGlycositeOverlaps(const std::vector<Residue*>& overlapResidues,
                                    const std::vector<Residue*>& glycositeResidues);
cds::Overlap countTotalOverlaps(const std::vector<std::vector<Residue*>>& overlapResidues,
                                const std::vector<std::vector<Residue*>>& glycositeResidues);
std::vector<size_t> determineSitesWithOverlap(int overlapTolerance,
                                              const std::vector<std::vector<Residue*>>& overlapResidues,
                                              const std::vector<std::vector<Residue*>>& glycositeResidues);
cds::Overlap wiggleSitesWithOverlaps(int overlapTolerance, int persistCycles, bool firstLinkageOnly, int interval,
                                     std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                     const std::vector<std::vector<Residue*>>& overlapResidues,
                                     const std::vector<std::vector<Residue*>>& glycositeResidues);
cds::Overlap randomDescent(bool monte_carlo, int persistCycles, int overlapTolerance,
                           std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                           const ResiduesByType& overlapResidues,
                           const std::vector<std::vector<Residue*>>& glycositeResidues);
bool dumbRandomWalk(int overlapTolerance, int maxCycles,
                    std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                    const std::vector<std::vector<Residue*>>& overlapResidues,
                    const std::vector<std::vector<Residue*>>& glycositeResidues);
void resolveOverlapsWithWiggler(bool randomize, int persistCycles, int overlapTolerance,
                                std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                const ResiduesByType& overlapResidues,
                                const std::vector<std::vector<Residue*>>& glycositeResidues);
#endif
