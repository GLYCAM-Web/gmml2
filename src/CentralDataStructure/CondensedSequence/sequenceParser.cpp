#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <sstream>
#include <string>
#include <cctype> // isdigit()
#include <vector>
#include <optional>

namespace
{
    using cds::ResidueType;
    using cdsCondensedSequence::AbstractSequence;
    using cdsCondensedSequence::ParsedResidueComponents;
    using cdsCondensedSequence::ResidueNode;

    const std::vector<char> openingBrackets {'(', '[', '{', '<'};
    const std::vector<char> closingBrackets {')', ']', '}', '>'};

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

    size_t matchingOpeningBracket(const std::string& str, size_t i)
    {
        std::vector<size_t> trace;
        char firstBracket = str[i];
        trace.reserve(64);
        size_t t = codeUtils::indexOf(closingBrackets, firstBracket);
        if (t > closingBrackets.size())
        {
            throw std::runtime_error("Not a closing bracket: '" + std::string {str[i]} + "' in: " + str);
        }
        trace.push_back(t);

        for (size_t n = i - 1; n < i; n--)
        {
            char c    = str[n];
            size_t ot = codeUtils::indexOf(openingBrackets, c);
            size_t ct = codeUtils::indexOf(closingBrackets, c);
            if (ct < closingBrackets.size())
            {
                trace.push_back(ct);
            }
            else if (ot < openingBrackets.size())
            {
                if (trace.empty())
                {
                    throw std::runtime_error("Could not parse brackets in: " + str);
                }
                else if (trace.back() == ot)
                {
                    if (trace.size() == 1)
                    {
                        return n;
                    }
                    trace.pop_back();
                }
            }
        }
        return str.length();
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
        auto mismatchingBracketError = [&](char first, char second, size_t position)
        {
            throw std::runtime_error("Did not find corresponding '" + std::string {first} + "' for character '" +
                                     std::string {second} + "' at position " + std::to_string(position + 1) +
                                     " of sequence: " + sequence);
        };
        for (size_t n = sequence.length() - 1; n < sequence.length(); n--)
        {
            size_t u = codeUtils::indexOf(openingBrackets, sequence[n]);
            if (u < openingBrackets.size())
            {
                mismatchingBracketError(closingBrackets[u], sequence[n], n);
            }
            size_t t = codeUtils::indexOf(closingBrackets, sequence[n]);
            if (t < closingBrackets.size())
            {
                size_t opening = matchingOpeningBracket(sequence, n);
                if (opening == sequence.length())
                {
                    mismatchingBracketError(openingBrackets[t], sequence[n], n);
                }
                n = opening;
            }
        }
        return true;
    }

    void parseLabelledInput(std::string inString)
    {
        throw std::runtime_error("Error: SequenceParser can't read labeled sequences yet: " + inString + "\n");
    }

    void updateResidue(AbstractSequence& data, size_t id, const ParsedResidueComponents& components)
    {
        data.nodes[id] = ResidueNode {components};
    }

    size_t addResidue(AbstractSequence& data, const ParsedResidueComponents& components)
    {
        data.nodes.emplace_back(ResidueNode {components});
        return graph::addNode(data.graph);
    }

    size_t addPlaceholderResidue(AbstractSequence& data)
    {
        return addResidue(data, {"", cds::ResidueType::Undefined, "", "", "", "", "", "", "", ""});
    }

    std::vector<size_t> recurveParseAlt(AbstractSequence& data, size_t parent, const std::string& sequence)
    {
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
                placeholder.reset();
                updateResidue(data, newRes, components);
                trail = "";
            }
            else
            {
                newRes = addResidue(data, components);
            }
            result.push_back(newRes);
            graph::addEdge(data.graph, {parent, newRes});
            return newRes;
        };
        size_t i         = sequence.length();
        size_t windowEnd = i;
        while (i > 0)
        {
            i--;
            if ((sequence[i] == '-') && ((windowEnd - i) > 5))
            { // dash and have read enough that there is a residue to save.
                parent    = save(substr(sequence, i + 2, windowEnd), parent);
                windowEnd = i + 2; // Get to the right side of the number e.g. 1-4
            }
            else if (sequence[i] == ']')
            { // Start of branch: recurve. Maybe save if have read enough.
                size_t bracketStart   = matchingOpeningBracket(sequence, i);
                std::string substring = substr(sequence, bracketStart + 1, i);
                if ((windowEnd - i) > 5)
                { // if not a derivative start and have read enough, save.
                    parent = save(substr(sequence, i + 1, windowEnd), parent);
                    recurveParseAlt(data, parent, substring);
                }
                else if ((windowEnd - i) > 1) // and not > 5
                {                             // Derivative, save trailing part of residue
                    trail       = substr(sequence, i + 1, windowEnd);
                    // reserve a residue to serve as a parent to derivatives, and update it later
                    placeholder = addPlaceholderResidue(data);
                    recurveParseAlt(data, placeholder.value(), substring);
                }
                else
                {
                    recurveParseAlt(data, parent, substring);
                }
                i         = bracketStart;
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
            else if (i == 0)
            { // Fin.
                save(substr(sequence, i, windowEnd), parent);
                return result;
            }
        }
        return result;
    }

    void parseCondensedSequence(AbstractSequence& data, const std::string& sequence)
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
        recurveParseAlt(data, terminal, substr(sequence, 0, i));
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
        size_t repeatCharacterEndLocation   = inputSequence.find_last_of('>');
        size_t repeatCharacterStartLocation = matchingOpeningBracket(inputSequence, repeatCharacterEndLocation);
        if (repeatCharacterStartLocation == inputSequence.length())
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
        size_t repeatStart         = matchingOpeningBracket(inputSequence, repeatEnd);
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
        if (before.find('>') != std::string::npos)
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

cdsCondensedSequence::AbstractSequence cdsCondensedSequence::parseSequence(std::string inputSequence)
{
    if (checkSequenceSanity(inputSequence))
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Sequence passed initial sanity checks for things like special characters or incorrect branching.\n");
    }
    AbstractSequence data;
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
        parseCondensedSequence(data, inputSequence);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "parseSequence complete");
    return data;
}
