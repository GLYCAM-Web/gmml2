#ifndef INCLUDE_PROGRAMS_GLYCOSYLATIONSITEFINDER_HPP
#define INCLUDE_PROGRAMS_GLYCOSYLATIONSITEFINDER_HPP

#include "include/fileType/pdb/pdbData.hpp"
#include "include/metadata/glycoprotein.hpp"

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
            const AminoAcidLinkTable& table, const pdb::PdbData& pdbData, const std::vector<size_t>& residues);
    } // namespace gpbuilder
} // namespace gmml

#endif