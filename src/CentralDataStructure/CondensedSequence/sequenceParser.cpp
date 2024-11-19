#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <sstream>
#include <string>
#include <cctype> // isdigit()
#include <vector>

using cdsCondensedSequence::ParsedResidue;

namespace
{
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

    ParsedResidue* saveResidue(std::vector<std::unique_ptr<ParsedResidue>>& savedResidues,
                               std::vector<std::string>& savedDerivatives, std::string residueString,
                               ParsedResidue* parent)
    {
        bool isDerivative = residueString.find('-') == std::string::npos;
        if (isDerivative)
        {
            savedDerivatives.push_back(residueString);
            return parent;
        }
        else
        {
            savedResidues.emplace_back(std::make_unique<ParsedResidue>(residueString, parent));
            ParsedResidue* newRes = savedResidues.back().get();
            for (auto& derivative : savedDerivatives)
            {
                savedResidues.emplace_back(std::make_unique<ParsedResidue>(derivative, newRes));
            }
            savedDerivatives.clear();
            return newRes;
        }
    }

    size_t recurveParseAlt(std::vector<std::unique_ptr<ParsedResidue>>& residues,
                           std::vector<std::string>& savedDerivatives, size_t i, const std::string& sequence,
                           ParsedResidue* parent)
    {
        auto save =
            [&residues, &savedDerivatives, &sequence](size_t windowStart, size_t windowEnd, ParsedResidue* parent)
        {
            std::string residueString = withoutBranches(substr(sequence, windowStart, windowEnd));
            return saveResidue(residues, savedDerivatives, residueString, parent);
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
                    i         = recurveParseAlt(residues, savedDerivatives, i, sequence, parent);
                    windowEnd = i;
                }
                else if ((windowEnd - i) > 1) // and not > 5
                {                             // Derivative. Note that windowEnd does not move.
                    i = recurveParseAlt(residues, savedDerivatives, i, sequence, parent);
                }
                else
                { // Return from branch and find new one. e.g. [Gal][Gal]
                    i         = recurveParseAlt(residues, savedDerivatives, i, sequence, parent);
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

    bool parseCondensedSequence(std::vector<std::unique_ptr<ParsedResidue>>& residues,
                                std::vector<std::string>& savedDerivatives, const std::string& sequence)
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
            residues.emplace_back(std::make_unique<ParsedResidue>(sequence.substr(i), cds::ResidueType::Sugar));
        }
        else
        { // e.g. DGlcpa1-OH
            residues.emplace_back(std::make_unique<ParsedResidue>(sequence.substr(i), cds::ResidueType::Aglycone));
        }
        ParsedResidue* terminal = residues.back().get();
        recurveParseAlt(residues, savedDerivatives, i, sequence, terminal);
        return true;
    }

    // Note Rob waved the wand and changed the format as of 2023-05-08. Changes below reflect that. Tails didn't change,
    // Heads yes and have become optional. Examples: DGlcpa1-2[4DGlcpa1-]<4>OH becomes
    // DGlcpa1-2DGlcpa1-4DGlcpa1-4DGlcpa1-4DGlcpa1-OH [4DGlcpa1-]<4>OH becomes DGlcpa1-4DGlcpa1-4DGlcpa1-4DGlcpa1-OH
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

void cdsCondensedSequence::parseSequence(cds::Molecule* molecule, std::string inputSequence)
{
    std::vector<std::string> savedDerivatives;
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
            std::vector<std::unique_ptr<ParsedResidue>> residues;
            parseCondensedSequence(residues, savedDerivatives, inputSequence);
            for (auto& res : residues)
            {
                molecule->addResidue(std::move(res));
            }
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "parseSequence complete");
    return;
}
