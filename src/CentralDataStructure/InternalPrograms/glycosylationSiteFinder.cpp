#include "includes/CentralDataStructure/InternalPrograms/glycosylationSiteFinder.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/MolecularMetadata/glycoprotein.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

using glycoproteinBuilder::GlycosylationSiteFinder;

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

    std::string GetSequenceContextAndDetermineTags(cds::Residue* residue, std::vector<std::string>& tags)
    {
        // Tags based on residue name only:
        std::string conjugationResidueCode     = glycoproteinMetadata::LookupCodeForAminoAcidName(residue->getName());
        std::vector<std::string> nameBasedTags = GetTagsForSequence(conjugationResidueCode, "", "");
        codeUtils::insertInto(tags, nameBasedTags);
        // Now look for context around residue via atom connectivity.
        std::string precedingContext         = "";
        cds::Residue* firstPrecedingNeighbor = cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(residue, "N");
        cds::Residue* secondPrecedingNeighbor = nullptr;
        if (firstPrecedingNeighbor)
        {
            precedingContext = glycoproteinMetadata::LookupCodeForAminoAcidName(firstPrecedingNeighbor->getName());
            secondPrecedingNeighbor =
                cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(firstPrecedingNeighbor, "N");
            if (secondPrecedingNeighbor)
            {
                precedingContext =
                    glycoproteinMetadata::LookupCodeForAminoAcidName(secondPrecedingNeighbor->getName()) +
                    precedingContext;
            }
        }
        std::string followingContext         = "";
        cds::Residue* firstFollowingNeighbor = cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(residue, "C");
        cds::Residue* secondFollowingNeighbor = nullptr;
        if (firstFollowingNeighbor)
        {
            std::string midResCode =
                glycoproteinMetadata::LookupCodeForAminoAcidName(firstFollowingNeighbor->getName());
            followingContext = midResCode;
            secondFollowingNeighbor =
                cdsSelections::FindNeighborResidueConnectedViaSpecificAtom(firstFollowingNeighbor, "C");
            if (secondFollowingNeighbor)
            {
                std::string lastResCode =
                    glycoproteinMetadata::LookupCodeForAminoAcidName(secondFollowingNeighbor->getName());
                followingContext += lastResCode;
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

GlycosylationSiteFinder::GlycosylationSiteFinder(std::vector<cds::Residue*> residues)
{
    for (auto& residue : residues)
    {
        std::string linkType = glycoproteinMetadata::LookupLinkTypeForAminoAcidName(residue->getName());
        if (!linkType.empty())
        {
            // std::cout << "Checking: " << residue->GetId() << "(" << residue->GetIndex() << ") with linkType: " <<
            // linkType << std::endl;
            std::vector<std::string> tags = {linkType};
            std::string context           = GetSequenceContextAndDetermineTags(residue, tags);
            pdb::ResidueId residueId      = residue->getId();
            table_.emplace_back(residueId.getChainId(), residueId.getNumber(), residueId.getInsertionCode(), context,
                                tags);
        }
    }
}

std::string GlycosylationSiteFinder::PrintTable()
{
    std::string output = "";
    for (auto& tableElement : table_)
    {
        output += tableElement.Print() + "\n";
    }
    return output;
}
