#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <sstream>
#include <string>
#include <cctype> // isdigit()
#include <vector>
#include <optional>

namespace
{
    using cds::ResidueType;
    using cdsCondensedSequence::BranchNode;
    using cdsCondensedSequence::NodeType;
    using cdsCondensedSequence::ParsedResidueComponents;
    using cdsCondensedSequence::ProbabilityNode;
    using cdsCondensedSequence::SequenceData;

    struct WeightedNodes
    {
        std::vector<std::vector<size_t>> chains;
        std::vector<double> weights;
        std::vector<bool> defaultWeight;
    };

    const std::string probabilityStr = "p=";

    std::string substr(const std::string& str, size_t windowStart, size_t windowEnd)
    {
        return str.substr(windowStart, (windowEnd - windowStart));
    }

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

    ParsedResidueComponents parseSugar(const std::string& fullString)
    { // E.g. DManpNAca1-4 . Isomer (D or L), residueName (ManNAc), ring type (f or p), configuration (a or b), linkage
        // (1-4)
        ParsedResidueComponents result;
        result.fullString         = fullString;
        std::string residueString = fullString;
        if (std::isdigit(fullString[0]))
        {
            residueString        = substr(fullString, 1, fullString.length() + 1);
            result.chainPosition = std::string {fullString[0]};
        }
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
    size_t seekBracketStart(const std::string& str, size_t i)
    {
        uint branches = 0;
        while (i > 0)
        {
            --i; // skip the initial ]
            if (str[i] == ']')
            {
                ++branches;
            }
            if (str[i] == '[')
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
        throw std::runtime_error("Could not find opening bracket in: " + str);
    }

    size_t firstOpeningBracket(const std::string& str, size_t pos)
    {
        for (size_t i = pos; i < str.size(); i++)
        {
            if (str[i] == '[')
            {
                return i;
            }
        }
        return str.length();
    }

    void parseLabelledInput(std::string inString)
    {
        throw std::runtime_error("Error: SequenceParser can't read labeled sequences yet: " + inString + "\n");
    }

    void updateResidue(SequenceData& data, size_t id, const ParsedResidueComponents& components)
    {
        data.nodes.fullString[id]        = components.fullString;
        data.nodes.type[id]              = components.type;
        data.nodes.name[id]              = components.name;
        data.nodes.linkage[id]           = components.linkage;
        data.nodes.chainPosition[id]     = components.chainPosition;
        data.nodes.ringType[id]          = components.ringType;
        data.nodes.configuration[id]     = components.configuration;
        data.nodes.isomer[id]            = components.isomer;
        data.nodes.preIsomerModifier[id] = components.preIsomerModifier;
        data.nodes.ringShape[id]         = components.ringShape;
        data.nodes.modifier[id]          = components.modifier;
    }

    size_t addNode(SequenceData& data, const NodeType& nodeType, const ParsedResidueComponents& components)
    {
        data.nodes.fullString.push_back(components.fullString);
        data.nodes.type.push_back(components.type);
        data.nodes.name.push_back(components.name);
        data.nodes.linkage.push_back(components.linkage);
        data.nodes.chainPosition.push_back(components.chainPosition);
        data.nodes.ringType.push_back(components.ringType);
        data.nodes.configuration.push_back(components.configuration);
        data.nodes.isomer.push_back(components.isomer);
        data.nodes.preIsomerModifier.push_back(components.preIsomerModifier);
        data.nodes.ringShape.push_back(components.ringShape);
        data.nodes.modifier.push_back(components.modifier);
        data.nodes.isInternal.push_back(false);
        data.nodes.isDerivative.push_back(
            codeUtils::contains({cds::ResidueType::Deoxy, cds::ResidueType::Derivative}, components.type));
        data.nodes.nodeType.push_back(nodeType);
        size_t id = graph::addNode(data.graph);
        return id;
    }

    size_t addResidue(SequenceData& data, const ParsedResidueComponents& components)
    {
        return addNode(data, cdsCondensedSequence::ResidueNode(), components);
    }

    size_t addPlaceholderResidue(SequenceData& data)
    {
        return addResidue(data, {});
    }

    void addLinkage(SequenceData& data, size_t child, size_t parent)
    {
        graph::addEdge(data.graph, {parent, child});
        data.edges.names.push_back("");
    }

    struct RecurveParseState
    {
        size_t parent;
        std::vector<size_t> chain;
    };

    RecurveParseState parseAny(SequenceData& data, size_t parent, const std::string& sequence);

    WeightedNodes parseWeightedNodes(SequenceData& data, size_t parent, const std::string& str)
    {
        std::vector<std::vector<size_t>> chains;
        std::vector<double> weights;
        std::vector<bool> defaultWeight;
        std::string errorStr = "Malformed syntax in: " + str;

        size_t i = str.length() - 1;

        bool pastFirst = false;
        while (i < str.length())
        {
            // the format is []or[]
            if (pastFirst && (str[i - 1] != 'o' || str[i] != 'r'))
            {
                throw std::runtime_error("'or' missing between options in: " + str);
            }
            if (pastFirst)
            {
                i -= 2;
            }
            pastFirst = true;
            if (str[i] != ']')
            {
                throw std::runtime_error(errorStr);
            }
            size_t localEnd         = i;
            size_t localStart       = seekBracketStart(str, localEnd);
            std::string localStr    = substr(str, localStart, localEnd + 1);
            size_t eqPosition       = localStr.length() - 1;
            double useDefaultWeight = false;
            for (size_t n = localStr.length() - 2; n > 0; n--)
            {
                char c = localStr[n];
                if (c == '=')
                {
                    eqPosition = n;
                    break;
                }
                else if (!(std::isdigit(c) || c == '.'))
                {
                    useDefaultWeight = true;
                    break;
                }
            }
            double weight = 1.0;
            if (!useDefaultWeight)
            {
                std::string weightStr        = substr(localStr, eqPosition + 1, localStr.length() - 1);
                std::optional<double> parsed = codeUtils::parseDouble(weightStr);
                if (!parsed.has_value())
                {
                    throw std::runtime_error(errorStr + " could not parse " + weightStr);
                }
                weight = parsed.value();
            }
            RecurveParseState state = parseAny(data, parent, substr(localStr, 1, eqPosition));
            weights.push_back(weight);
            defaultWeight.push_back(useDefaultWeight);
            chains.push_back(state.chain);
            i = localStart - 1;
        }

        return {chains, weights, defaultWeight};
    }

    RecurveParseState parseProbability(SequenceData& data, const size_t initialParent, const std::string& sequence)
    {
        size_t openingBracket             = firstOpeningBracket(sequence, 1);
        std::string numberStr             = substr(sequence, probabilityStr.length(), openingBracket);
        std::optional<double> probability = codeUtils::parseDouble(numberStr);
        if (!probability.has_value())
        {
            throw std::runtime_error("Couldn't parse probability '" + numberStr + "' in: " + sequence);
        }
        size_t nodeId = addNode(data, ProbabilityNode {}, {sequence});
        addLinkage(data, nodeId, initialParent);
        WeightedNodes weighted = parseWeightedNodes(data, nodeId, substr(sequence, openingBracket, sequence.length()));
        // we want at least one option
        if (weighted.chains.empty())
        {
            throw std::runtime_error("Probability lacks options in: " + sequence);
        }
        // if there are multiple, we want them all to have explicit proportions
        else if (weighted.chains.size() > 1)
        {
            for (bool def : weighted.defaultWeight)
            {
                if (def)
                {
                    throw std::runtime_error("Multiple probability options must contain relative proportions in: " +
                                             sequence);
                }
            }
        }
        bool isDerivative               = data.nodes.isDerivative[weighted.chains[0][0]];
        data.nodes.isDerivative[nodeId] = isDerivative;
        for (std::vector<size_t>& chain : weighted.chains)
        {
            if (data.nodes.isDerivative[chain[0]] != isDerivative)
            {
                throw std::runtime_error("Probability options must all be sugars or all be derivatives in: " +
                                         sequence);
            }
        }
        ProbabilityNode& node = std::get<ProbabilityNode>(data.nodes.nodeType[nodeId]);
        node.probability      = probability.value();
        if (probability < 0.0 || probability > 1.0)
        {
            throw std::runtime_error("Probability value must be between 0 and 1 in: " + sequence);
        }
        node.heads.reserve(weighted.chains.size());
        node.tails.reserve(weighted.chains.size());
        for (auto& chain : weighted.chains)
        {
            node.heads.push_back(chain.front());
            node.tails.push_back(chain.back());
        }
        node.weights = weighted.weights;
        return {nodeId, {nodeId}};
    }

    RecurveParseState parseBrackets(SequenceData& data, const size_t initialParent, const std::string& sequence)
    {
        std::string str = substr(sequence, 1, sequence.length() - 1);
        if (codeUtils::startsWith(str, probabilityStr))
        {
            return parseProbability(data, initialParent, str);
        }
        else if (!str.empty() && isdigit(str.back()))
        { // branch, digit at end is linkage
            std::string linkage = {str.back()};
            size_t branchId     = addNode(data, BranchNode {linkage, 0}, {sequence});
            addLinkage(data, branchId, initialParent);
            // parse without linkage number
            RecurveParseState state = parseAny(data, branchId, substr(str, 0, str.length() - 1));
            BranchNode& node        = std::get<BranchNode>(data.nodes.nodeType[branchId]);
            node.head               = state.chain[0];
            return {initialParent, {}};
        }
        else
        {
            return parseAny(data, initialParent, str);
        }
    }

    RecurveParseState parseList(SequenceData& data, const size_t initialParent, const std::string& sequence)
    {
        size_t parent = initialParent;
        std::vector<size_t> result;
        result.reserve(32);
        std::optional<size_t> placeholder = {};
        std::string trail                 = "";
        auto save                         = [&](const std::string& str, size_t parent)
        {
            ParsedResidueComponents components =
                parseResidueStringIntoComponents(str + trail, cds::ResidueType::Undefined);
            size_t newRes = -1;
            if (placeholder.has_value())
            { // we've reserved a residue spot to update
                newRes = placeholder.value();
                updateResidue(data, newRes, components);
                placeholder.reset();
                trail = "";
            }
            else
            {
                newRes = addResidue(data, components);
            }
            result.push_back(newRes);
            addLinkage(data, newRes, parent);
            return newRes;
        };
        size_t i         = sequence.length();
        size_t windowEnd = i;
        while (i > 0)
        {
            i--;
            if (i == 0)
            { // Fin.
                save(substr(sequence, i, windowEnd), parent);
                return {initialParent, result};
            }
            else if (std::isdigit(sequence[i]) && (sequence[i - 1] == '-') && ((windowEnd - i) + trail.length() > 6))
            { // dash and have read enough that there is a residue to save.
                parent    = save(substr(sequence, i, windowEnd), parent);
                windowEnd = i;
            }
            else if (sequence[i] == ']')
            { // Start of branch: recurve. Maybe save if have read enough.
                size_t startBracket   = seekBracketStart(sequence, i);
                std::string substring = substr(sequence, startBracket, i + 1);
                if ((windowEnd - i) > 5)
                { // if not a derivative start and have read enough, save.
                    parent                  = save(substr(sequence, i + 1, windowEnd), parent);
                    RecurveParseState state = parseBrackets(data, parent, substring);
                    codeUtils::insertInto(result, state.chain);
                    parent = state.parent;
                }
                else if ((windowEnd - i) > 1) // and not > 5
                {                             // Derivative, save trailing part of residue
                    trail       = substr(sequence, i + 1, windowEnd);
                    // reserve a residue to serve as a parent to derivatives, and update it later
                    placeholder = addPlaceholderResidue(data);
                    parseBrackets(data, placeholder.value(), substring);
                }
                else
                {
                    RecurveParseState state = parseBrackets(data, parent, substring);
                    codeUtils::insertInto(result, state.chain);
                    parent = state.parent;
                }
                i         = startBracket;
                windowEnd = i;
            }
            else if (sequence[i] == '[')
            {
                throw std::runtime_error("Unclosed bracket in: " + sequence);
            }
            else if (sequence[i] == ',')
            { // , in derivative list
                save(substr(sequence, i + 1, windowEnd), parent);
                windowEnd = i;
            }
        }
        return {initialParent, result};
    }

    RecurveParseState parseAny(SequenceData& data, size_t parent, const std::string& sequence)
    {
        if (!sequence.empty() && sequence.back() == ']' && seekBracketStart(sequence, sequence.length() - 1) == 0)
        {
            return parseBrackets(data, parent, sequence);
        }
        else
        {
            return parseList(data, parent, sequence);
        }
    }

    void parseCondensedSequence(SequenceData& data, const std::string& sequence)
    {
        // Reading from the rightmost end of the string, get the aglycone first.
        size_t i = (sequence.find_last_of('-') + 1);
        if (substr(sequence, sequence.length() - 2, sequence.length()) == "OH")
        {
            addResidue(data, parseAglycone(sequence.substr(sequence.length() - 2)));
            parseAny(data, 0, substr(sequence, 0, sequence.length() - 2));
        }
        else if (isdigit(sequence[i])) // Indicates ano-ano and not e.g. Sugar-OME or -ROH etc
        {                              // e.g. DGlcpa1-2DFrufb
            ++i;                       // ano-ano
            if (sequence[i] == ']')    // e.g. DGlcpa1-2[LFucpa1-1]DFrufb also ano-ano, but with branching
            {
                ++i;
            }
            addResidue(data, parseResidueStringIntoComponents(sequence.substr(i), cds::ResidueType::Sugar));
            size_t terminal = data.graph.nodes.size() - 1;
            parseAny(data, terminal, substr(sequence, 0, i));
        }
        else
        {
            addResidue(data, parseAglycone(sequence.substr(i)));
            parseAny(data, 0, substr(sequence, 0, i));
        }
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
        size_t repeatStart         = seekBracketStart(inputSequence, repeatEnd);
        std::string before         = inputSequence.substr(0, repeatStart);
        std::string repeat         = substr(inputSequence, repeatStart + 1, repeatEnd);
        std::string after          = inputSequence.substr(repeatCharacterEndLocation + 1);
        std::string newInputString = before;
        // Check if using the old nomenclature. i.e. DGlcpa1-[4DGlcpa1-]<3>
        bool oldNomenclature       = !before.empty() && before.back() == '-';
        bool numberFirstEndOpen    = !repeat.empty() && std::isdigit(repeat[0]) && repeat.back() == '-';
        if (numberFirstEndOpen && !oldNomenclature)
        { // New nomenclature, i.e DGlcpa1-4[4DGlcpa1-]<3>
            // firstRepeat does not have the e.g. 4 in 4DGlcpa1-
            newInputString += repeat.substr(1, repeat.length());
            for (int n = 0; n < numberRepeats - 1; n++)
            {
                newInputString += repeat;
            }
        }
        else
        { // all other cases, repeat as-is
            for (int n = 0; n < numberRepeats; n++)
            {
                newInputString += repeat;
            }
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
