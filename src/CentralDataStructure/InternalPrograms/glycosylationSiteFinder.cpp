#include "includes/CentralDataStructure/InternalPrograms/glycosylationSiteFinder.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/MolecularMetadata/glycoprotein.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>
#include <vector>

namespace
{
    std::vector<std::string> GetTagsForSequence(const std::string& firstRes, const std::string& midRes,
                                                const std::string& lastRes)
    {
        std::vector<std::string> foundTags;
        // N-links
        if ((firstRes == "N") && (midRes != "P"))
        {
            if ((lastRes == "T") || (lastRes == "S"))
            {
                foundTags.push_back("nLinkLikely");
            }
            if (lastRes == "C")
            {
                foundTags.push_back("nLinkCysSequon");
            }
        }
        // O-links // The midRes == "" is funky, but if I'm coming in here with more than just firstRes, I already have
        // this tag from the initial call.
        if (((firstRes == "S") || (firstRes == "T")) && (midRes == ""))
        {
            foundTags.push_back("oLinkLikely");
        }
        if ((firstRes == "Y") && (midRes == ""))
        {
            foundTags.push_back("oLinkTyr");
        }
        return foundTags;
    }

    std::string GetSequenceContextAndDetermineTags(const pdb::PdbData& pdbData, cds::Residue* residue,
                                                   std::vector<std::string>& tags)
    {
        const glycoproteinMetadata::AminoAcidLinkTable& table = glycoproteinMetadata::defaultAminoAcidLinkTable();
        auto getCode                                          = [&table](const std::string& name)
        {
            size_t index = codeUtils::indexOf(table.names, name);
            if (index == table.names.size())
            {
                throw std::runtime_error("Error: no amino acid table entry found for: " + name);
            }
            return table.codes[index];
        };
        // Tags based on residue name only:
        std::string conjugationResidueCode     = getCode(residue->getName());
        std::vector<std::string> nameBasedTags = GetTagsForSequence(conjugationResidueCode, "", "");
        codeUtils::insertInto(tags, nameBasedTags);
        // Now look for context around residue via atom connectivity.
        std::string precedingContext = "";
        cds::Residue* firstPrecedingNeighbor =
            cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(pdbData, residue, "N");
        cds::Residue* secondPrecedingNeighbor = nullptr;
        if (firstPrecedingNeighbor)
        {
            precedingContext = getCode(firstPrecedingNeighbor->getName());
            secondPrecedingNeighbor =
                cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(pdbData, firstPrecedingNeighbor, "N");
            if (secondPrecedingNeighbor)
            {
                precedingContext = getCode(secondPrecedingNeighbor->getName()) + precedingContext;
            }
        }
        std::string followingContext = "";
        cds::Residue* firstFollowingNeighbor =
            cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(pdbData, residue, "C");
        cds::Residue* secondFollowingNeighbor = nullptr;
        if (firstFollowingNeighbor)
        {
            std::string midResCode = getCode(firstFollowingNeighbor->getName());
            followingContext       = midResCode;
            secondFollowingNeighbor =
                cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(pdbData, firstFollowingNeighbor, "C");
            if (secondFollowingNeighbor)
            {
                std::string lastResCode = getCode(secondFollowingNeighbor->getName());
                followingContext        += lastResCode;
                std::vector<std::string> contextTags =
                    GetTagsForSequence(conjugationResidueCode, midResCode, lastResCode);
                codeUtils::insertInto(tags, contextTags);
            }
        }
        std::string context = precedingContext + "_" + conjugationResidueCode + "_" + followingContext;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Context: " + context);
        return context;
    }
} // namespace

namespace glycoproteinBuilder
{
    std::vector<GlycosylationSiteInfo> createGlycosylationSiteTable(const pdb::PdbData& pdbData,
                                                                    const std::vector<cds::Residue*>& residues)
    {
        const glycoproteinMetadata::AminoAcidLinkTable& table = glycoproteinMetadata::defaultAminoAcidLinkTable();
        std::vector<GlycosylationSiteInfo> result;
        for (auto& residue : residues)
        {
            size_t index = codeUtils::indexOf(table.names, residue->getName());
            if (index < table.names.size() && table.linkTypes[index] != "")
            {
                size_t residueId              = codeUtils::indexOf(pdbData.objects.residues, residue);
                std::vector<std::string> tags = {table.linkTypes[index]};
                std::string context           = GetSequenceContextAndDetermineTags(pdbData, residue, tags);
                result.push_back(GlycosylationSiteInfo {pdbData.residues.chainIds[residueId],
                                                        std::to_string(residue->getNumber()),
                                                        pdbData.residues.insertionCodes[residueId], context, tags});
            }
        }
        return result;
    }
} // namespace glycoproteinBuilder
