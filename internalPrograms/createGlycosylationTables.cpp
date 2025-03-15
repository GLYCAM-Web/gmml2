#include "includes/CentralDataStructure/InternalPrograms/glycosylationSiteFinder.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp"
#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/structuredFiles.hpp"
#include "includes/version.h"

#include <string>
#include <vector>
#include <functional>
#include <iostream>

int main(int argc, char* argv[])
{
    enum ARGUMENTS
    {
        FILENAME,
        HELP,
        VERSION,
        FORMAT
    };

    enum OUTPUT_FORMAT
    {
        LIST,
        CSV,
        TXT,
        HTML
    };

    using codeUtils::ArgReq;
    using codeUtils::ArgType;
    // keeping the same order between the enum and strings lets us assign format by index in the vector
    std::vector<std::string> knownFormats {"list", "csv", "txt", "html"};
    std::vector<codeUtils::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed, FILENAME, "", ' ', "filename"},
        {ArgReq::optional, ArgType::flag, HELP, "help", 'h', ""},
        {ArgReq::optional, ArgType::flag, VERSION, "version", 'v', ""},
        {ArgReq::optional, ArgType::option, FORMAT, "format", 'f', codeUtils::join("|", knownFormats)},
    };
    std::string programName = codeUtils::programName(argv);

    codeUtils::Arguments arguments;
    try
    {
        arguments = codeUtils::readArguments(argc, argv, argumentDefinitions);
        if (codeUtils::contains<int>(arguments.ids, HELP))
        {
            std::cout << codeUtils::helpString(programName, argumentDefinitions);
            std::cout << "\n"
                      << "For more information, see https://github.com/GLYCAM-Web/gmml2\n";
            std::exit(0);
        }
        if (codeUtils::contains<int>(arguments.ids, VERSION))
        {
            std::cout << "Glycoprotein Builder Table & GMML version " << GMML_VERSION << "\n";
            std::exit(0);
        }
        codeUtils::validateArguments(arguments, argumentDefinitions);
    }
    catch (const std::runtime_error& error)
    {
        std::cout << "error in program arguments\n";
        std::cout << error.what() << "\n";
        std::cout << "\n";
        std::cout << codeUtils::helpString(programName, argumentDefinitions);
        std::exit(1);
    }

    std::string inputFileName = "";
    OUTPUT_FORMAT format      = LIST;
    for (const auto& arg : arguments.args)
    {
        switch (arg.id)
        {
            case ARGUMENTS::FILENAME:
                {
                    inputFileName = arg.value;
                    break;
                }
            case ARGUMENTS::FORMAT:
                {
                    size_t index = codeUtils::indexOf(knownFormats, arg.value);
                    if (index >= knownFormats.size())
                    {
                        throw std::runtime_error("Unknown format: '" + arg.value + "', expected one of " +
                                                 codeUtils::join(", ", knownFormats));
                    }
                    format = OUTPUT_FORMAT(index);
                    break;
                }
            default:
                break;
        }
    }
    using glycoproteinBuilder::GlycosylationSiteInfo;
    pdb::PdbFile inputFile(inputFileName);
    std::vector<cds::Residue*> residues = pdb::getResidues(inputFile.getAssemblies());
    cds::bondAtomsAndResiduesByDistance(residues);
    std::vector<GlycosylationSiteInfo> table = glycoproteinBuilder::createGlycosylationSiteTable(residues);
    std::vector<std::string> header {"Chain", "ResidueNumber", "InsertionCode", "SequenceContext", "Tags"};

    std::function<std::vector<std::string>(const GlycosylationSiteInfo&)> toRow = [](const GlycosylationSiteInfo& info)
    {
        return std::vector<std::string> {info.chain, info.residueNumber, info.insertionCode, info.sequenceContext,
                                         codeUtils::join(" ", info.tags)};
    };
    std::vector<std::vector<std::string>> rows = codeUtils::vectorMap(toRow, table);
    codeUtils::TextTable textTable {header, rows};
    switch (format)
    {
        case LIST:
            for (auto& a : table)
            {
                std::cout << "Chain: " << a.chain << "\nResidueNumber: " << a.residueNumber
                          << "\nInsertionCode: " + a.insertionCode << "\nSequenceContext: " << a.sequenceContext
                          << "\nTags:\n"
                          << codeUtils::join("\n", a.tags) << "\n\n";
            }
            break;
        case CSV:
            codeUtils::toCsv(std::cout, ",", textTable);
            break;
        case TXT:
            codeUtils::toTxt(std::cout, {textTable});
            break;
        case HTML:
            codeUtils::toHtml(std::cout, {textTable});
            break;
    }
    return 0;
}
