#ifndef INCLUDES_INTERNALPROGRAMS_GLYCOSYLATIONSITEFINDER_HPP
#define INCLUDES_INTERNALPROGRAMS_GLYCOSYLATIONSITEFINDER_HPP

#include "include/readers/Pdb/pdbData.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
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
    } // namespace gpbuilder
} // namespace gmml

#endif