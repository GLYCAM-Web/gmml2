#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <sstream>
#include <string>
#include <cctype> // isdigit()
#include <vector>

namespace
{
    using cds::ResidueType;
    using cdsCondensedSequence::ParsedResidueComponents;
    using cdsCondensedSequence::SequenceData;

    std::pair<std::string, std::string> ringShapeAndModifier(const std::string& modifier)
    { // E.g. LIdopA(4C1)a1-4 with modifier "A(4C1)", which here gets broken into ring shape "4C1" and modifier "A".
        size_t leftParenthesisPosition  = modifier.find('(');
        size_t rightParenthesisPosition = modifier.find(')');

        if ((leftParenthesisPosition == std::string::npos) || (rightParenthesisPosition == std::string::npos))
        { // If there isn't a ring shape declared.
            return {"", modifier};
        }
        else
        { // Assumes it's always at end of modifiers
            return {modifier.substr(leftParenthesisPosition), modifier.substr(0, leftParenthesisPosition)};
        }
    }

    ParsedResidueComponents parseAglycone(const std::string& residueString)
    {
        return {residueString, ResidueType::Aglycone, residueString, "", "", "", "", "", "", ""};
    }

    ParsedResidueComponents parseDerivative(const std::string& residueString)
    { // A derivative e.g. 3S, 6Me. Linkage followed by residue name. No configuration.
        std::string name      = residueString.substr(1); // From position 1 to the end.
        cds::ResidueType type = (name == "D") || (name == "H") ? ResidueType::Deoxy : ResidueType::Derivative;
        std::string linkage   = residueString.substr(0, 1);
        return {residueString, type, name, linkage, "", "", "", "", "", ""};
    }

    ParsedResidueComponents parseSugar(const std::string& residueString)
    { // E.g. DManpNAca1-4 . Isomer (D or L), residueName (ManNAc), ring type (f or p), configuration (a or b), linkage
        // (1-4)
        ParsedResidueComponents result;
        result.fullString    = residueString;
        result.type          = ResidueType::Sugar;
        // Assumptions
        size_t isomerStart   = 0; // e.g. DGal, DGlc, LIdo
        size_t residueStart  = 0; // e.g. Gal, Glc, Ido
        size_t modifierStart = 3; // e.g. NAc, A, A(1C4)
        // Checks
        if ((residueString.find("LDmanpHep") != std::string::npos) ||
            (residueString.find("DDmanpHep") != std::string::npos)) // Exceptions for these weirdos
        {
            result.preIsomerModifier = residueString.substr(0, 1);
            ++residueStart;
            ++modifierStart;
            ++isomerStart;
        }
        std::string isomer = residueString.substr(isomerStart, 1);
        if ((isomer == "D") || (isomer == "L"))
        {
            result.isomer = isomer;
            ++residueStart;  // 1 if normal e.g. DGlc..
            ++modifierStart; // 4 if normal.
        }
        result.name          = residueString.substr(residueStart, 3);
        size_t ringPosition  = (residueStart + 3);
        std::string ringType = residueString.substr(ringPosition, 1);
        if ((ringType == "p") || (ringType == "f"))
        {
            result.ringType = ringType;
            ++modifierStart; // 5 if normal.
        }
        // Find the dash, read around it.
        size_t dashPosition = residueString.find('-');
        if (dashPosition == std::string::npos) // There is no -. e.g. Fru in DGlcpa1-2DFrufb
        {
            dashPosition = residueString.size() + 1;
        }
        else
        {
            result.linkage = residueString.substr((dashPosition - 1), 3);
        }
        size_t modifierEnd = dashPosition - 2; // They are 2 apart if no modifier i.e. DGlcpa1-2, the "a1" size is 2
        std::string configuration = residueString.substr(dashPosition - 2, 1);
        if ((configuration == "a" || configuration == "b"))
        {
            result.configuration = configuration;
        }
        else
        {
            modifierEnd++; // e.g. if ano is missing DGlcpA1-OH
        }
        // Find any special modifiers e.g. NAc, Gc, A in DGlcpAa1-OH NAc in DGlcpNAca1-2
        size_t modifierLength = (modifierEnd - modifierStart);
        if (modifierLength > 100)
        {
            std::string message = "This is a non-standard residue string that gmml can't handle: " + residueString;
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        if (modifierLength > 0 && modifierLength < 100)
        {
            std::string modifier                           = residueString.substr(modifierStart, modifierLength);
            std::pair<std::string, std::string> ringAndMod = ringShapeAndModifier(modifier);
            result.ringShape                               = ringAndMod.first;
            result.modifier                                = ringAndMod.second;
        }
        return result;
    }

    ParsedResidueComponents parseResidueStringIntoComponents(const std::string& residueString,
                                                             ResidueType specifiedType)
    {
        if ((residueString.find('-') != std::string::npos) || (specifiedType == ResidueType::Sugar))
        {
            return parseSugar(residueString);
        }
        else if (isdigit(residueString[0]))
        {
            return parseDerivative(residueString);
        }
        else
        { // Dunno.
            std::string message = "Error: we can't parse this residue: \"" + residueString + "\"";
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
    }

    std::string substr(const std::string& str, size_t windowStart, size_t windowEnd)
    {
        return str.substr(windowStart, (windowEnd - windowStart));
    }

    bool checkSequenceSanity(const std::string& sequence)
    {
        if (sequence.empty())
        {
            throw std::runtime_error("Error: sequence is empty:>>>" + sequence + "<<<");
        }
        if (sequence.find("cake") != std::string::npos)
        {
            throw std::runtime_error("Error: the cake is a lie:>>>" + sequence + "<<<");
        }
        if (sequence.find(" ") != std::string::npos)
        {
            throw std::runtime_error("Error: sequence contains a space:>>>" + sequence + "<<<");
        }
        std::vector<char> badChars = {'\'', '_', '+', '"', '`'};
        for (auto& badChar : badChars)
        {
            if (sequence.find(badChar) != std::string::npos)
            {
                std::string s(1, badChar); // convert to string
                throw std::runtime_error("Error: sequence cannot contain this:\'" + s + "\':>>>" + sequence + "<<<");
            }
        }
        size_t a = std::count(sequence.begin(), sequence.end(), '[');
        size_t b = std::count(sequence.begin(), sequence.end(), ']');
        if (a != b)
        {
            throw std::runtime_error("Error: the number of [ doesn't match the number of ]. Bad branch in :>>>" +
                                     sequence + "<<<");
        }
        return true;
    }

    // There may be branches which also use [], so need to check for those and find the [ that starts the repeat
    unsigned int seekRepeatStart(const std::string& inputSequence, unsigned int i)
    {
        int branches = 0;

        while (i > 0)
        {
            --i; // skip the initial ]
            if (inputSequence[i] == ']')
            {
                ++branches;
            }
            if (inputSequence[i] == '[')
            {
                if (branches == 0)
                {
                    return i;
                }
                else
                {
                    --branches;
                }
            }
        }
        throw std::runtime_error("Did not find corresponding '[' in repeat unit of repeating sequence: " +
                                 inputSequence);
    }

    void parseLabelledInput(std::string inString)
    {
        throw std::runtime_error("Error: SequenceParser can't read labeled sequences yet: " + inString + "\n");
    }

    std::string withoutBranches(const std::string& residueString)
    {
        //  Splice out anything within [ and ].
        size_t branch_start = residueString.find_first_of('[');
        if (branch_start != std::string::npos)
        {
            size_t branch_finish  = residueString.find_last_of(']') + 1;
            std::string firstPart = residueString.substr(0, branch_start);
            std::string lastPart  = residueString.substr(branch_finish);
            return firstPart + lastPart;
        }
        else
        {
            return residueString;
        }
    }

    size_t addResidue(SequenceData& data, const ParsedResidueComponents& components)
    {
        data.residues.fullString.push_back(components.fullString);
        data.residues.type.push_back(components.type);
        data.residues.name.push_back(components.name);
        data.residues.linkage.push_back(components.linkage);
        data.residues.ringType.push_back(components.ringType);
        data.residues.configuration.push_back(components.configuration);
        data.residues.isomer.push_back(components.isomer);
        data.residues.preIsomerModifier.push_back(components.preIsomerModifier);
        data.residues.ringShape.push_back(components.ringShape);
        data.residues.modifier.push_back(components.modifier);
        data.residues.isInternal.push_back(false);
        size_t id = graph::addNode(data.graph);
        return id;
    }

    void addLinkage(SequenceData& data, size_t resA, size_t resB)
    {
        bool isSugar         = data.residues.type[resA] == cds::ResidueType::Sugar;
        bool isChildDeoxy    = data.residues.type[resA] == cds::ResidueType::Deoxy;
        std::string edgeName = (isSugar ? data.residues.configuration[resA] : "") + data.residues.linkage[resA];
        graph::addEdge(data.graph, {resB, resA});
        // It remains internal if it's already been made internal, or if the child is a deoxy.
        data.residues.isInternal[resB] = (!isChildDeoxy || data.residues.isInternal[resB]);
        data.edges.names.push_back(edgeName);
    }

    size_t saveResidue(SequenceData& data, std::vector<size_t>& savedDerivatives, std::string residueString,
                       size_t parent)
    {
        bool isDerivative = residueString.find('-') == std::string::npos;
        if (isDerivative)
        {
            size_t derivative =
                addResidue(data, parseResidueStringIntoComponents(residueString, cds::ResidueType::Undefined));
            savedDerivatives.push_back(derivative);
            return parent;
        }
        else
        {
            size_t newRes =
                addResidue(data, parseResidueStringIntoComponents(residueString, cds::ResidueType::Undefined));
            addLinkage(data, newRes, parent);
            for (size_t derivative : savedDerivatives)
            {
                addLinkage(data, derivative, newRes);
            }
            savedDerivatives.clear();
            return newRes;
        }
    }

    size_t recurveParseAlt(SequenceData& data, std::vector<size_t>& savedDerivatives, size_t i,
                           const std::string& sequence, size_t parent)
    {
        auto save = [&](size_t windowStart, size_t windowEnd, size_t parent)
        {
            std::string residueString = withoutBranches(substr(sequence, windowStart, windowEnd));
            return saveResidue(data, savedDerivatives, residueString, parent);
        };
        size_t windowEnd = i;
        while (i > 0)
        {
            i--;
            if ((sequence[i] == '-') && ((windowEnd - i) > 5))
            { // dash and have read enough that there is a residue to save.
                parent    = save(i + 2, windowEnd, parent);
                windowEnd = i + 2; // Get to the right side of the number e.g. 1-4
            }
            else if (sequence[i] == ']')
            { // Start of branch: recurve. Maybe save if have read enough.
                if ((windowEnd - i) > 5)
                { // if not a derivative start and have read enough, save.
                    parent    = save(i + 1, windowEnd, parent);
                    i         = recurveParseAlt(data, savedDerivatives, i, sequence, parent);
                    windowEnd = i;
                }
                else if ((windowEnd - i) > 1) // and not > 5
                {                             // Derivative. Note that windowEnd does not move.
                    i = recurveParseAlt(data, savedDerivatives, i, sequence, parent);
                }
                else
                { // Return from branch and find new one. e.g. [Gal][Gal]
                    i         = recurveParseAlt(data, savedDerivatives, i, sequence, parent);
                    windowEnd = i;
                }
            }
            else if (sequence[i] == '[')
            { // End of branch
                save(i + 1, windowEnd, parent);
                return i;
            }
            else if (sequence[i] == ',')
            { // , in derivative list
                parent    = save(i + 1, windowEnd, parent);
                windowEnd = i;
            }
            else if (i == 0)
            { // Fin.
                parent = save(i, windowEnd, parent);
            }
        }
        return i;
    }

    void parseCondensedSequence(SequenceData& data, const std::string& sequence)
    {
        // Reading from the rightmost end of the string, get the aglycone first.
        size_t i = (sequence.find_last_of('-') + 1);
        if (isdigit(sequence[i]))   // Indicates ano-ano and not e.g. Sugar-OME or -ROH etc
        {                           // e.g. DGlcpa1-2DFrufb
            ++i;                    // ano-ano
            if (sequence[i] == ']') // e.g. DGlcpa1-2[LFucpa1-1]DFrufb also ano-ano, but with branching
            {
                ++i;
            }
            addResidue(data, parseResidueStringIntoComponents(sequence.substr(i), cds::ResidueType::Sugar));
        }
        else
        { // e.g. DGlcpa1-OH
            addResidue(data, parseAglycone(sequence.substr(i)));
        }
        size_t terminal = data.graph.nodes.size() - 1;
        std::vector<size_t> savedDerivatives;
        recurveParseAlt(data, savedDerivatives, i, sequence, terminal);
    }

    // Note Rob waved the wand and changed the format as of 2023-05-08. Changes below reflect that. Tails didn't change,
    // Heads yes and have become optional. Examples:
    // DGlcpa1-2[4DGlcpa1-]<4>OH becomes DGlcpa1-2DGlcpa1-4DGlcpa1-4DGlcpa1-4DGlcpa1-OH
    // [4DGlcpa1-]<4>OH becomes DGlcpa1-4DGlcpa1-4DGlcpa1-4DGlcpa1-OH
    // DGlcpa1-3[4DGlcpa1-]<9>2DManpa1-[4DGalpNAca1-]<4>OH // Multiple repeats
    // DGlcpa1-2[4DGlcpa1-3DManpa1-]<9>OH // Disacc repeats
    // DGlcpa1-2[4DGlcpa1-3[DAllpb1-2]DManpa1-]<9>OH // Branched repeats, unavailable on legacy
    // Repeats within repeats?
    // DGlcpa1-4[4[DGalpa1-3]DGlcpa1-3[DAllpb1-2]DManpa1-]<3>OH Repeats with branches on the leftmost Residue?
    // DGlcpa1-[4DGlcpa1-3DManpa1-]<3>4DGalpa1-OH // Tails that aren't OH
    std::string parseRepeatingUnits(const std::string& inputSequence)
    {
        size_t repeatCharacterEndLocation   = inputSequence.find_first_of('>');
        size_t repeatCharacterStartLocation = inputSequence.find_first_of('<');
        if (repeatCharacterStartLocation == std::string::npos)
        {
            throw std::runtime_error("No '<' found in sequence with repeating syntax symbol '>' : " + inputSequence);
        }
        // Get the number of repeats:
        int numberRepeats = 0;
        try
        {
            size_t numberStart       = repeatCharacterStartLocation + 1;
            std::string stringNumber = substr(inputSequence, numberStart, repeatCharacterEndLocation);
            numberRepeats            = std::stoi(stringNumber);
        }
        catch (...) // if e.g. stoi throws
        {
            throw std::runtime_error("Number of repeating units not specified correctly in repeating unit: " +
                                     inputSequence);
        }
        size_t repeatEnd = repeatCharacterStartLocation - 1;
        // Ensure next char is the ] of the repeating unit
        if (inputSequence[repeatEnd] != ']')
        {
            throw std::runtime_error("Missing or incorrect usage of ']' in repeating sequence: " + inputSequence);
        }
        // Ok now go find the position of the start of the repeating unit, considering branches
        size_t repeatStart       = seekRepeatStart(inputSequence, repeatEnd);
        std::string beforeRepeat = inputSequence.substr(0, repeatStart);
        // Check if using the old nomenclature. i.e. DGlcpa1-[4DGlcpa1-]<4> was old, DGlcpa1-4[4DGlcpa1-]<3> is new.
        if (inputSequence.find("-[") != std::string::npos)
        {
            std::string message = "Old repeat syntax detected as \"-[\" found in " + inputSequence;
            gmml::log(__LINE__, __FILE__, gmml::INF, message);
            //  I think I'll just copy the e.g. 4 to make it look like the new syntax.
            std::string numberToInsert = inputSequence.substr(repeatStart + 1, 1);
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Adding this before repeat unit to make it be new syntax: " + numberToInsert);
            beforeRepeat += numberToInsert;
        }

        std::string firstRepeat = inputSequence.substr(
            repeatStart + 2, (repeatEnd - repeatStart - 2)); // firstRepeat does not have the e.g. 4 in 4DGlcpa1-
        std::string repeat = substr(inputSequence, repeatStart + 1, repeatEnd);
        std::string after  = inputSequence.substr(repeatCharacterEndLocation + 1);
        std::string newInputString;
        newInputString += beforeRepeat;
        newInputString += firstRepeat;
        numberRepeats--; // firstRepeat added already.
        for (int j = 1; j <= numberRepeats; ++j)
        {
            newInputString += repeat;
        }
        newInputString += after;
        // Check if there are more repeating units and deal with them recursively:
        if (after.find('>') != std::string::npos)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Sequence with some repeats processed: " + newInputString);
            return parseRepeatingUnits(newInputString);
        }
        else
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Sequence with all repeats processed: " + newInputString);
            return newInputString;
        }
    }
} // namespace

std::vector<size_t> cdsCondensedSequence::edgesSortedByLink(const SequenceData& sequence,
                                                            const std::vector<size_t>& edgeIds)
{
    auto residueLink = [&](size_t n)
    {
        return cdsCondensedSequence::getLink(sequence.residues.type[n], sequence.residues.linkage[n]);
    };
    std::function<bool(const size_t&, const size_t&)> compare = [&](const size_t& n, const size_t& k)
    {
        return residueLink(sequence.graph.edgeNodes[n][1]) > residueLink(sequence.graph.edgeNodes[k][1]);
    };

    return codeUtils::sortedBy(compare, edgeIds);
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::reordered(const SequenceData& sequence)
{
    size_t residueCount           = sequence.residues.name.size();
    std::vector<size_t> edgeOrder = edgesSortedByLink(sequence, codeUtils::indexVector(sequence.graph.edges));
    std::vector<std::array<size_t, 2>> reorderedEdges = codeUtils::indicesToValues(sequence.graph.edgeNodes, edgeOrder);

    size_t current                   = 0;
    std::vector<size_t> residueOrder = {current};
    residueOrder.reserve(residueCount);
    std::vector<bool> traversed = codeUtils::indicesToBools(residueCount, residueOrder);
    auto fromNode               = [&](size_t node)
    {
        std::vector<size_t> result;
        result.reserve(reorderedEdges.size());
        for (auto& edge : reorderedEdges)
        {
            if (edge[0] == node)
            {
                result.push_back(edge[1]);
            }
        }
        return codeUtils::reverse(result);
    };
    std::vector<size_t> nodesToTraverse;
    codeUtils::insertInto(nodesToTraverse, fromNode(current));
    while (!nodesToTraverse.empty())
    {
        current = nodesToTraverse.back();
        nodesToTraverse.pop_back();
        if (!traversed[current])
        {
            residueOrder.push_back(current);
            codeUtils::insertInto(nodesToTraverse, fromNode(current));
            traversed[current] = true;
        }
    }

    std::vector<size_t> missed = codeUtils::boolsToIndices(codeUtils::vectorNot(traversed));
    if (missed.size() > 0)
    {
        throw std::runtime_error("Error: sequence graph not fully connected");
    }

    std::vector<size_t> invertedResidueOrder(residueCount, -1);
    for (size_t n = 0; n < residueOrder.size(); n++)
    {
        invertedResidueOrder[residueOrder[n]] = n;
    }

    graph::Database resultGraph;

    for (size_t n = 0; n < residueOrder.size(); n++)
    {
        graph::addNode(resultGraph);
    }

    std::vector<size_t> edgeIndexOrder = codeUtils::indexVector(sequence.graph.edges);

    for (size_t n : edgeIndexOrder)
    {
        const std::array<size_t, 2>& edge = sequence.graph.edgeNodes[n];
        graph::addEdge(resultGraph, {invertedResidueOrder[edge[0]], invertedResidueOrder[edge[1]]});
    }

    ResidueData residues = {codeUtils::indicesToValues(sequence.residues.fullString, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.type, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.name, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.linkage, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringType, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.configuration, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isomer, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.preIsomerModifier, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringShape, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.modifier, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isInternal, residueOrder)};

    return SequenceData {resultGraph, residues, sequence.edges};
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::parseSequence(std::string inputSequence)
{
    SequenceData data;
    if (inputSequence.find('<') != std::string::npos)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Found repeating unit in input\n");
        inputSequence = parseRepeatingUnits(inputSequence);
    }
    if (inputSequence.find(';') != std::string::npos)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Found labels in input\n");
        parseLabelledInput(inputSequence);
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Parsing unlabelled input sequence:\n" + inputSequence + "\n");
        if (checkSequenceSanity(inputSequence))
        {
            gmml::log(
                __LINE__, __FILE__, gmml::INF,
                "Sequence passed initial sanity checks for things like special characters or incorrect branching.\n");
            parseCondensedSequence(data, inputSequence);
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "parseSequence complete");
    return data;
}
