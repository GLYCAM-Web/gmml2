#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOSYLATIONSITEFINDER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOSYLATIONSITEFINDER_HPP

#include <string>
#include <vector>
#include "includes/CentralDataStructure/residue.hpp"

namespace glycoproteinBuilder
{
    struct GlycosylationSiteInfo
    {
        std::string chain;
        std::string residueNumber;
        std::string insertionCode;
        std::string sequenceContext;
        std::vector<std::string> tags; // e.g. oLink, nLink, sequon, cysteineSequon, all, etc
    };

    std::vector<GlycosylationSiteInfo> createGlycosylationSiteTable(const std::vector<cds::Residue*>& residues);
} // namespace glycoproteinBuilder

#endif