#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/bondByDistance.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/metadata/elements.hpp"
#include "include/metadata/glycoprotein.hpp"
#include "include/programs/glycosylationSiteFinder.hpp"
#include "include/util/arguments.hpp"
#include "include/util/containers.hpp"
#include "include/util/strings.hpp"
#include "include/util/structuredFiles.hpp"
#include "include/version.h"

#include <functional>
#include <iostream>
#include <string>
#include <vector>

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
        TXT,
        CSV,
        HTML
    };

    using namespace gmml;
    using util::ArgReq;
    using util::ArgType;
    // keeping the same order between the enum and strings lets us assign format by index in the vector
    std::vector<std::string> knownFormats {"txt", "csv", "html"};
    std::vector<util::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed, FILENAME, "", ' ', "filename"},
        {ArgReq::optional, ArgType::flag, HELP, "help", 'h', ""},
        {ArgReq::optional, ArgType::flag, VERSION, "version", 'v', ""},
        {ArgReq::optional, ArgType::option, FORMAT, "format", 'f', util::join("|", knownFormats)},
    };
    std::string programName = util::programName(argv);

    util::Arguments arguments;
    try
    {
        arguments = util::readArguments(argc, argv, argumentDefinitions);
        if (util::contains<int>(arguments.ids, HELP))
        {
            std::cout << util::helpString(programName, argumentDefinitions);
            std::cout << "\n"
                      << "For more information, see https://github.com/GLYCAM-Web/gmml2\n";
            std::exit(0);
        }
        if (util::contains<int>(arguments.ids, VERSION))
        {
            std::cout << "Glycoprotein Builder Table & GMML2 version " << GMML_VERSION << "\n";
            std::exit(0);
        }
        util::validateArguments(arguments, argumentDefinitions);
    }
    catch (const std::runtime_error& error)
    {
        std::cout << "error in program arguments\n";
        std::cout << error.what() << "\n";
        std::cout << "\n";
        std::cout << util::helpString(programName, argumentDefinitions);
        std::exit(1);
    }

    std::string inputFileName = "";
    OUTPUT_FORMAT format = TXT;
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
                    size_t index = util::indexOf(knownFormats, arg.value);
                    if (index >= knownFormats.size())
                    {
                        throw std::runtime_error(
                            "Unknown format: '" + arg.value + "', expected one of " + util::join(", ", knownFormats));
                    }
                    format = OUTPUT_FORMAT(index);
                    break;
                }
            default:
                break;
        }
    }
    pdb::PdbFile pdbFile = pdb::toPdbFile(inputFileName, {pdb::InputType::modelsAsMolecules, false});
    pdb::PdbData& pdbData = pdbFile.data;
    util::SparseVector<double> elementRadii = vanDerWaalsRadii();
    const std::vector<Sphere> atomBounds =
        assembly::toAtomBounds(elementRadii, pdbData.atoms.elements, pdbData.atoms.coordinates);
    const assembly::Graph graph = assembly::createAssemblyGraph(pdbData.assembly.indices, pdbData.assembly.atomGraph);
    const assembly::Bounds bounds = assembly::toAssemblyBounds(graph, atomBounds);
    pdb::bondAtomsAndResiduesByDistance(pdbData, bounds);
    const AminoAcidLinkTable linkTable = defaultAminoAcidLinkTable();
    std::vector<gpbuilder::GlycosylationSiteInfo> table =
        gpbuilder::createGlycosylationSiteTable(linkTable, pdbData, util::indexVector(residueCount(pdbData.assembly)));
    std::vector<std::string> header {"Chain", "ResidueNumber", "InsertionCode", "SequenceContext", "Tags"};

    std::function<std::vector<std::string>(const gpbuilder::GlycosylationSiteInfo&)> toRow =
        [](const gpbuilder::GlycosylationSiteInfo& info)
    {
        return std::vector<std::string> {
            info.chain, info.residueNumber, info.insertionCode, info.sequenceContext, util::join(" ", info.tags)};
    };
    std::vector<std::vector<std::string>> rows = util::vectorMap(toRow, table);
    util::TextTable textTable {header, rows};
    switch (format)
    {
        case TXT:
            util::toTxt(std::cout, {textTable});
            break;
        case CSV:
            util::toCsv(std::cout, ",", textTable);
            break;
        case HTML:
            util::toHtml(std::cout, {textTable});
            break;
    }
    return 0;
}
