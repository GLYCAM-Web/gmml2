#include "include/programs/glycosylationSiteFinder.hpp"

#include "include/CentralDataStructure/Selections/residueSelections.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/metadata/glycoprotein.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        namespace
        {
            std::vector<std::string> GetTagsForSequence(
                const std::string& firstRes, const std::string& midRes, const std::string& lastRes)
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
                // O-links // The midRes == "" is funky, but if I'm coming in here with more than just firstRes, I
                // already have this tag from the initial call.
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

            size_t firstNeighborInOtherResidue(
                const pdb::PdbData& pdbData, const graph::Graph& atomGraph, size_t atomId)
            {
                size_t residueId = pdbData.indices.atomResidue[atomId];
                for (size_t n : atomGraph.nodes.nodeAdjacencies[atomId])
                {
                    if (pdbData.indices.atomResidue[n] != residueId)
                    {
                        return n;
                    }
                }
                return pdbData.indices.atomCount;
            }

            // Not having Atom know which Residue it is in makes this funky. Make a decision about whether that happens
            // or not.
            size_t findNeighborResidueConnectedViaSpecificAtom(
                const pdb::PdbData& pdbData,
                const graph::Graph& atomGraph,
                size_t queryResidueId,
                const std::string& queryAtomName)
            {
                size_t atomCount = pdbData.indices.atomCount;
                size_t residueCount = pdbData.indices.residueCount;
                std::string queryResidueStr = pdb::residueStringId(pdbData, queryResidueId);
                size_t queryAtom = pdb::findResidueAtom(pdbData, queryResidueId, queryAtomName);
                if (queryAtom == atomCount)
                {
                    util::log(
                        __LINE__,
                        __FILE__,
                        util::WAR,
                        "Warning: An atom named " + queryAtomName + " not found in residue: " + queryResidueStr);
                    return residueCount;
                }
                size_t foreignAtomNeighbor = firstNeighborInOtherResidue(pdbData, atomGraph, queryAtom);
                if (foreignAtomNeighbor == atomCount)
                {
                    util::log(
                        __LINE__,
                        __FILE__,
                        util::WAR,
                        "Warning: Did not find foreign neighbors for an atom named " + queryAtomName +
                            " in residue: " + queryResidueStr);
                    return residueCount;
                }
                return pdbData.indices.atomResidue[foreignAtomNeighbor];
            }

            std::string GetSequenceContextAndDetermineTags(
                const AminoAcidLinkTable& table,
                const pdb::PdbData& pdbData,
                const graph::Graph& atomGraph,
                size_t residueId,
                std::vector<std::string>& tags)
            {
                size_t residueCount = pdbData.indices.residueCount;
                auto getCode = [&](size_t residueId)
                {
                    const std::string& name = pdbData.residues.names[residueId];
                    size_t index = util::indexOf(table.names, name);
                    if (index == table.names.size())
                    {
                        throw std::runtime_error("Error: no amino acid table entry found for: " + name);
                    }
                    return table.codes[index];
                };
                // Tags based on residue name only:
                std::string conjugationResidueCode = getCode(residueId);
                std::vector<std::string> nameBasedTags = GetTagsForSequence(conjugationResidueCode, "", "");
                util::insertInto(tags, nameBasedTags);
                // Now look for context around residue via atom connectivity.
                std::string precedingContext = "";
                size_t firstPrecedingNeighbor =
                    findNeighborResidueConnectedViaSpecificAtom(pdbData, atomGraph, residueId, "N");
                if (firstPrecedingNeighbor < residueCount)
                {
                    precedingContext = getCode(firstPrecedingNeighbor);
                    size_t secondPrecedingNeighbor =
                        findNeighborResidueConnectedViaSpecificAtom(pdbData, atomGraph, firstPrecedingNeighbor, "N");
                    if (secondPrecedingNeighbor < residueCount)
                    {
                        precedingContext = getCode(secondPrecedingNeighbor) + precedingContext;
                    }
                }
                std::string followingContext = "";
                size_t firstFollowingNeighbor =
                    findNeighborResidueConnectedViaSpecificAtom(pdbData, atomGraph, residueId, "C");
                if (firstFollowingNeighbor < residueCount)
                {
                    std::string midResCode = getCode(firstFollowingNeighbor);
                    followingContext = midResCode;
                    size_t secondFollowingNeighbor =
                        findNeighborResidueConnectedViaSpecificAtom(pdbData, atomGraph, firstFollowingNeighbor, "C");
                    if (secondFollowingNeighbor < residueCount)
                    {
                        std::string lastResCode = getCode(secondFollowingNeighbor);
                        followingContext += lastResCode;
                        std::vector<std::string> contextTags =
                            GetTagsForSequence(conjugationResidueCode, midResCode, lastResCode);
                        util::insertInto(tags, contextTags);
                    }
                }
                std::string context = precedingContext + "_" + conjugationResidueCode + "_" + followingContext;
                util::log(__LINE__, __FILE__, util::INF, "Context: " + context);
                return context;
            }
        } // namespace

        std::vector<GlycosylationSiteInfo> createGlycosylationSiteTable(
            const pdb::PdbData& pdbData, const std::vector<size_t>& residues)
        {
            const AminoAcidLinkTable& table = defaultAminoAcidLinkTable();
            graph::Graph atomGraph = graph::identity(pdbData.atomGraph);
            std::vector<GlycosylationSiteInfo> result;
            for (size_t residueId : residues)
            {
                size_t index = util::indexOf(table.names, pdbData.residues.names[residueId]);
                if (index < table.names.size() && table.linkTypes[index] != "")
                {
                    std::vector<std::string> tags = {table.linkTypes[index]};
                    std::string context =
                        GetSequenceContextAndDetermineTags(table, pdbData, atomGraph, residueId, tags);
                    result.push_back(GlycosylationSiteInfo {
                        pdbData.residues.chainIds[residueId],
                        std::to_string(pdbData.residues.numbers[residueId]),
                        pdbData.residues.insertionCodes[residueId],
                        context,
                        tags});
                }
            }
            return result;
        }
    } // namespace gpbuilder
} // namespace gmml
