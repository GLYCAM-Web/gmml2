#include "include/sequence/sequenceParser.hpp"

#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/sequence/sequenceTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/parsing.hpp"

#include <cctype> // isdigit()
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace gmml
{
    namespace sequence
    {
        namespace
        {
            const std::vector<char> openingBrackets {'(', '[', '{', '<'};
            const std::vector<char> closingBrackets {')', ']', '}', '>'};

            const size_t noDerivatives = 0;

            std::pair<std::string, std::string> ringShapeAndModifier(const std::string& modifier)
            { // E.g. LIdopA(4C1)a1-4 with modifier "A(4C1)", which here gets broken into ring shape "4C1" and modifier
              // "A".
                size_t leftParenthesisPosition = modifier.find('(');
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
                return {residueString, ResidueType::Aglycone, residueString, "", "", "", "", "", "", "", ""};
            }

            ParsedResidueComponents parseDerivative(const std::string& residueString)
            { // A derivative e.g. 3S, 6Me. Linkage followed by residue name. No configuration.
                std::string name = residueString.substr(1); // From position 1 to the end.
                ResidueType type = (name == "D") || (name == "H") ? ResidueType::Deoxy : ResidueType::Derivative;
                std::string linkage = residueString.substr(0, 1);
                return {residueString, type, name, linkage, "", "", "", "", "", "", ""};
            }

            ParsedResidueComponents parseSugar(std::string residueString)
            { // E.g. DManpNAca1-4 . Isomer (D or L), residueName (ManNAc), ring type (f or p), configuration (a or b),
              // linkage
                // (1-4)
                ParsedResidueComponents result;
                result.fullString = residueString;
                result.type = ResidueType::Sugar;
                if (std::isdigit(residueString[0]))
                {
                    result.defaultHeadPosition = residueString.substr(0, 1);
                    residueString = residueString.substr(1);
                }
                // Assumptions
                size_t isomerStart = 0;   // e.g. DGal, DGlc, LIdo
                size_t residueStart = 0;  // e.g. Gal, Glc, Ido
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
                result.name = residueString.substr(residueStart, 3);
                size_t ringPosition = (residueStart + 3);
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
                size_t modifierEnd =
                    dashPosition - 2; // They are 2 apart if no modifier i.e. DGlcpa1-2, the "a1" size is 2
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
                    std::string message =
                        "This is a non-standard residue string that gmml can't handle: " + residueString;
                    util::log(__LINE__, __FILE__, util::ERR, message);
                    throw std::runtime_error(message);
                }
                if (modifierLength > 0 && modifierLength < 100)
                {
                    std::string modifier = residueString.substr(modifierStart, modifierLength);
                    std::pair<std::string, std::string> ringAndMod = ringShapeAndModifier(modifier);
                    result.ringShape = ringAndMod.first;
                    result.modifier = ringAndMod.second;
                }
                return result;
            }

            bool isAglycone(const std::string& str) { return util::contains({"OH", "OME", "OtBu"}, str); }

            SequenceNode parseResidueString(const std::string& residueString, size_t derivativeList)
            {
                if (isAglycone(residueString))
                {
                    return AglyconeNode {parseAglycone(residueString)};
                }
                else if (isdigit(residueString[0]) && residueString.length() <= 3)
                {
                    return DerivativeNode {parseDerivative(residueString)};
                }
                else if (residueString.length() >= 5)
                {
                    return MonosaccharideNode {parseSugar(residueString), derivativeList};
                }
                else
                {
                    throw std::runtime_error("Could not parse this residue: " + residueString);
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
                size_t t = util::indexOf(closingBrackets, firstBracket);
                if (t > closingBrackets.size())
                {
                    throw std::runtime_error("Not a closing bracket: '" + std::string {str[i]} + "' in: " + str);
                }
                trace.push_back(t);

                for (size_t n = i - 1; n < i; n--)
                {
                    char c = str[n];
                    size_t ot = util::indexOf(openingBrackets, c);
                    size_t ct = util::indexOf(closingBrackets, c);
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
                        throw std::runtime_error(
                            "Error: sequence cannot contain this:\'" + s + "\':>>>" + sequence + "<<<");
                    }
                }
                auto mismatchingBracketError = [&](char first, char second, size_t position)
                {
                    throw std::runtime_error(
                        "Did not find corresponding '" + std::string {first} + "' for character '" +
                        std::string {second} + "' at position " + std::to_string(position + 1) +
                        " of sequence: " + sequence);
                };
                for (size_t n = sequence.length() - 1; n < sequence.length(); n--)
                {
                    size_t u = util::indexOf(openingBrackets, sequence[n]);
                    if (u < openingBrackets.size())
                    {
                        mismatchingBracketError(closingBrackets[u], sequence[n], n);
                    }
                    size_t t = util::indexOf(closingBrackets, sequence[n]);
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

            size_t addNode(AbstractSequence& data, const std::string& str, const SequenceNode& node)
            {
                size_t id = data.nodes.size();
                data.nodes.emplace_back(node);
                data.fullString.push_back(str);
                return id;
            }

            std::vector<size_t> parseChain(AbstractSequence& data, const std::string& sequence)
            {
                std::vector<size_t> result;
                result.reserve(32);
                size_t derivativeList = noDerivatives;
                std::string trail = "";
                auto save = [&](const std::string& str)
                {
                    std::string fullStr = str + trail;
                    SequenceNode node = parseResidueString(fullStr, derivativeList);
                    trail = "";
                    derivativeList = noDerivatives;
                    size_t newRes = addNode(data, fullStr, node);
                    result.push_back(newRes);
                };
                size_t i = sequence.length();
                size_t windowEnd = i;
                while (i > 0)
                {
                    i--;
                    std::string after = substr(sequence, i + 1, windowEnd);
                    if ((sequence[i] == '-') && isAglycone(after))
                    {
                        save(after);
                        windowEnd = i + 1;
                    }
                    else if ((sequence[i] == '-') && (windowEnd - i > 5))
                    {
                        save(substr(sequence, i + 2, windowEnd));
                        windowEnd = i + 2; // Get to the right side of the number e.g. 1-4
                    }
                    else if (util::contains({']', '>'}, sequence[i]))
                    { // Start of branch: recurve. Maybe save if have read enough.
                        bool isSquare = sequence[i] == ']';
                        bool isAngular = sequence[i] == '>';
                        size_t bracketStart = matchingOpeningBracket(sequence, i);
                        std::string substringWithBrackets = substr(sequence, bracketStart, i + 1);
                        std::string substring = substr(sequence, bracketStart + 1, i);
                        int unsaved = (windowEnd - i - 1);
                        if (isSquare && unsaved > 0 && unsaved <= 4 && !isAglycone(after))
                        { // Derivative, save trailing part of residue
                            trail = after;
                            std::vector<size_t> derivatives = parseChain(data, substring);
                            derivativeList = addNode(data, substringWithBrackets, DerivativeListNode {derivatives});
                        }
                        else
                        {
                            if (unsaved > 0)
                            {
                                save(substr(sequence, i + 1, windowEnd));
                            }
                            if (isAngular)
                            {
                                // Note Rob waved the wand and changed the format as of 2023-05-08. Changes below
                                // reflect that. Tails didn't change, Heads yes and have become optional. Examples:
                                // DGlcpa1-2[4DGlcpa1-]<4>OH becomes DGlcpa1-2DGlcpa1-4DGlcpa1-4DGlcpa1-4DGlcpa1-OH
                                // [4DGlcpa1-]<4>OH becomes DGlcpa1-4DGlcpa1-4DGlcpa1-4DGlcpa1-OH
                                // DGlcpa1-3[4DGlcpa1-]<9>2DManpa1-[4DGalpNAca1-]<4>OH // Multiple repeats
                                // DGlcpa1-2[4DGlcpa1-3DManpa1-]<9>OH // Disacc repeats
                                // DGlcpa1-2[4DGlcpa1-3[DAllpb1-2]DManpa1-]<9>OH // Branched repeats, unavailable on
                                // legacy Repeats within repeats?
                                // DGlcpa1-4[4[DGalpa1-3]DGlcpa1-3[DAllpb1-2]DManpa1-]<3>OH Repeats with branches on the
                                // leftmost Residue? DGlcpa1-[4DGlcpa1-3DManpa1-]<3>4DGalpa1-OH // Tails that aren't OH
                                std::optional<uint> repeats = util::parseUint(substring);
                                if (!repeats.has_value() || repeats.value() == 0)
                                {
                                    throw std::runtime_error(
                                        "Number of repeating units not specified correctly in repeating unit: " +
                                        substringWithBrackets);
                                }
                                if (bracketStart == 0 || sequence[bracketStart - 1] != ']')
                                {
                                    throw std::runtime_error(
                                        "Repeat marker " + substringWithBrackets +
                                        " must be directly preceded by ']' in: " + sequence);
                                }
                                i = bracketStart - 1;
                                bracketStart = matchingOpeningBracket(sequence, i);
                                std::string repeatPart = substringWithBrackets;
                                std::string substringWithBrackets = substr(sequence, bracketStart, i + 1);
                                std::string substring = substr(sequence, bracketStart + 1, i);
                                std::vector<size_t> constituents = parseChain(data, substring);
                                size_t chain = addNode(data, substring, ChainNode {constituents});
                                size_t repeat = addNode(
                                    data, substringWithBrackets + repeatPart, RepeatNode {repeats.value(), chain});
                                result.push_back(repeat);
                            }
                            else
                            {
                                char last = substring[substring.length() - 1];
                                std::string position = "";
                                if (std::isdigit(last))
                                { // branch position is explicit
                                    position = std::string {last};
                                    substring = substring.substr(0, substring.length() - 1);
                                }
                                std::vector<size_t> constituents = parseChain(data, substring);
                                size_t chain = addNode(data, substring, ChainNode {constituents});
                                size_t branch = addNode(data, substringWithBrackets, BranchNode {position, chain});
                                result.push_back(branch);
                            }
                        }
                        i = bracketStart;
                        windowEnd = i;
                    }
                    else if (sequence[i] == '[')
                    {
                        throw std::runtime_error("Unclosed bracket in: " + sequence);
                    }
                    else if (sequence[i] == ',')
                    { // , in derivative list
                        save(substr(sequence, i + 1, windowEnd));
                        windowEnd = i;
                    }
                    else if (i == 0)
                    { // Fin.
                        save(substr(sequence, i, windowEnd));
                        return result;
                    }
                }
                return result;
            }

            void parseCondensedSequence(AbstractSequence& data, const std::string& sequence)
            {
                // all sugars will contain a list of derivatives,
                // thus we let the first node be an empty list for convenience
                addNode(data, "", DerivativeListNode {{}});
                std::vector<size_t> constituents = parseChain(data, substr(sequence, 0, sequence.length()));
                data.root = addNode(data, sequence, ChainNode {constituents});
            }
        } // namespace

        AbstractSequence parseSequence(std::string inputSequence)
        {
            if (checkSequenceSanity(inputSequence))
            {
                util::log(
                    __LINE__,
                    __FILE__,
                    util::INF,
                    "Sequence passed initial sanity checks for things like special characters or incorrect "
                    "branching.\n");
            }
            AbstractSequence data;
            if (inputSequence.find(';') != std::string::npos)
            {
                util::log(__LINE__, __FILE__, util::INF, "Found labels in input\n");
                throw std::runtime_error(
                    "Error: SequenceParser can't read labeled sequences yet: " + inputSequence + "\n");
            }
            else
            {
                util::log(__LINE__, __FILE__, util::INF, "Parsing unlabelled input sequence:\n" + inputSequence + "\n");
                parseCondensedSequence(data, inputSequence);
            }
            util::log(__LINE__, __FILE__, util::INF, "parseSequence complete");
            return data;
        }
    } // namespace sequence
} // namespace gmml
