#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOSYLATIONSITEFINDER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOSYLATIONSITEFINDER_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"

#include <string>
#include <vector>

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

    std::vector<GlycosylationSiteInfo> createGlycosylationSiteTable(
        const pdb::PdbData& pdbData, const std::vector<size_t>& residues);
} // namespace glycoproteinBuilder

#endif