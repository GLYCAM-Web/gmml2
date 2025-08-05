#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/fileType/pdb/bondByDistance.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/arguments.hpp"
#include "include/util/containers.hpp"
#include "include/util/parsing.hpp"
#include "include/util/strings.hpp"
#include "include/util/structuredFiles.hpp"
#include "include/version.h"

#include <cmath>
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
        FORMAT,
        MODE,
        ORDER,
        OVERLAPTOLERANCE
    };

    enum OUTPUT_FORMAT
    {
        CSV,
        TXT,
        HTML
    };

    enum OUTPUT_MODE
    {
        MOLECULES,
        RESIDUES,
        ATOMS,
        CONTACTS
    };

    enum OUTPUT_ORDER
    {
        INDEX,
        OVERLAP,
        POTENTIAL
    };

    using namespace gmml;
    using util::ArgReq;
    using util::ArgType;
    // keeping the same order between the enum and strings lets us assign format by index in the vector
    std::vector<std::string> knownFormats {"csv", "txt", "html"};
    std::vector<std::string> modes {"molecules", "residues", "atoms", "contacts"};
    std::vector<std::string> orders {"index", "overlap", "potential"};
    std::vector<util::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed, FILENAME, "", ' ', "filename"},
        {ArgReq::optional, ArgType::flag, HELP, "help", 'h', ""},
        {ArgReq::optional, ArgType::flag, VERSION, "version", 'v', ""},
        {ArgReq::optional, ArgType::option, FORMAT, "format", 'f', util::join("|", knownFormats)},
        {ArgReq::optional, ArgType::option, MODE, "mode", 'm', util::join("|", modes)},
        {ArgReq::optional, ArgType::option, ORDER, "order", 'o', util::join("|", orders)},
        {ArgReq::optional, ArgType::option, OVERLAPTOLERANCE, "tolerance", 't', "value"}
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
            std::cout << "PDB Analyzer & GMML2 version " << GMML_VERSION << "\n";
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
    OUTPUT_MODE mode = CONTACTS;
    OUTPUT_ORDER order = INDEX;
    double overlapTolerance = 0.0;
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
            case ARGUMENTS::MODE:
                {
                    size_t index = util::indexOf(modes, arg.value);
                    if (index >= modes.size())
                    {
                        throw std::runtime_error(
                            "Unknown mode: '" + arg.value + "', expected one of " + util::join(", ", modes));
                    }
                    mode = OUTPUT_MODE(index);
                    break;
                }
            case ARGUMENTS::ORDER:
                {
                    size_t index = util::indexOf(orders, arg.value);
                    if (index >= orders.size())
                    {
                        throw std::runtime_error(
                            "Unknown order: '" + arg.value + "', expected one of " + util::join(", ", orders));
                    }
                    order = OUTPUT_ORDER(index);
                    break;
                }
            case ARGUMENTS::OVERLAPTOLERANCE:
                {
                    std::optional<double> opt = util::parseDouble(arg.value);
                    if (opt.has_value())
                    {
                        overlapTolerance = opt.value();
                    }
                    else
                    {
                        throw std::runtime_error("Not a double: " + arg.value);
                    }
                    break;
                }
            default:
                break;
        }
    }
    pdb::PdbFile pdbFile = pdb::toPdbFile(inputFileName, pdb::modelsAsMolecules);
    pdb::PdbData& data = pdbFile.data;
    util::SparseVector<double> elementRadii = vanDerWaalsRadii();
    const PotentialTable potential = potentialTable(elementRadii, foundElements(data.atoms.elements));
    const std::vector<Sphere> atomBounds =
        assembly::toAtomBounds(elementRadii, data.atoms.elements, data.atoms.coordinates);
    const assembly::Graph graph = assembly::createAssemblyGraph(data.indices, data.atomGraph);
    const assembly::Bounds bounds = assembly::toAssemblyBounds(graph, atomBounds);

    pdb::bondAtomsAndResiduesByDistance(data, bounds);
    auto residueAtomsCloseToEdge = assembly::atomsCloseToResidueEdges(graph);

    std::vector<std::string> residueTypeNames(ResidueTypeCount, "");
    residueTypeNames[ResidueType::Aglycone] = "aglycone";
    residueTypeNames[ResidueType::Deoxy] = "deoxy";
    residueTypeNames[ResidueType::Derivative] = "derivative";
    residueTypeNames[ResidueType::Protein] = "protein";
    residueTypeNames[ResidueType::ProteinCappingGroup] = "protein cap";
    residueTypeNames[ResidueType::Solvent] = "solvent";
    residueTypeNames[ResidueType::Sugar] = "sugar";
    residueTypeNames[ResidueType::Undefined] = "unknown";

    std::function<std::vector<std::string>(const size_t&)> moleculeRow = [&](size_t moleculeId)
    {
        std::string chain = data.molecules.chainIds[moleculeId];
        std::vector<bool> includedTypes(ResidueTypeCount, false);
        for (ResidueType type : util::indicesToValues(data.residues.types, moleculeResidues(data.indices, moleculeId)))
        {
            includedTypes[type] = true;
        }
        std::string types = util::join(", ", util::boolsToValues(residueTypeNames, includedTypes));
        std::string residueCount = std::to_string(assembly::moleculeResidues(data.indices, moleculeId).size());
        std::string atomCount = std::to_string(assembly::moleculeAtoms(data.indices, moleculeId).size());
        return std::vector<std::string> {chain, residueCount, atomCount, types};
    };
    std::vector<std::string> moleculeHeader {"chain", "residues", "atoms", "residue types"};

    std::function<std::vector<std::string>(const size_t&)> residueRow = [&](size_t residueId)
    {
        std::string chain = data.residues.chainIds[residueId];
        std::string name = data.residues.names[residueId];
        std::string number = std::to_string(data.residues.numbers[residueId]);
        std::string type = residueTypeNames[data.residues.types[residueId]];
        std::string atomCount = std::to_string(assembly::residueAtoms(data.indices, residueId).size());
        return std::vector<std::string> {chain, name, number, type, atomCount};
    };
    std::vector<std::string> residueHeader {"chain", "name", "number", "type", "atom count"};

    std::function<std::vector<std::string>(const size_t&)> atomRow = [&](size_t atomId)
    {
        size_t residueId = data.indices.atomResidue[atomId];
        std::string chain = data.residues.chainIds[residueId];
        std::string residueName = data.residues.names[residueId];
        std::string residueNumber = std::to_string(data.residues.numbers[residueId]);
        std::string name = data.atoms.names[atomId];
        std::string number = std::to_string(data.atoms.numbers[atomId]);
        std::string element = elementName(data.atoms.elements[atomId]);
        return std::vector<std::string> {chain, residueName, residueNumber, name, number, element};
    };
    std::vector<std::string> atomHeader {"chain", "residue name", "residue number", "name", "number", "element"};

    struct AtomContact
    {
        size_t atomA;
        size_t atomB;
        double overlap;
        double potential;
    };

    std::function<std::vector<AtomContact>(size_t, size_t)> contactRows = [&](size_t residueA, size_t residueB)
    {
        std::vector<AtomContact> result;
        bool bothProtein = (data.residues.types[residueA] == ResidueType::Protein) &&
                           (data.residues.types[residueB] == ResidueType::Protein);
        bool eitherWater =
            util::contains({data.residues.names[residueA], data.residues.names[residueB]}, std::string("HOH"));
        if (!bothProtein && !eitherWater &&
            spheresOverlap(overlapTolerance, bounds.residues[residueA], bounds.residues[residueB]))
        {
            const std::vector<size_t>& atomsA = graph.residues.nodes.constituents[residueA];
            const std::vector<size_t>& atomsB = graph.residues.nodes.constituents[residueB];
            std::vector<bool> nonIgnoredA(atomsA.size(), true);
            std::vector<bool> nonIgnoredB(atomsB.size(), true);
            size_t adjacency = util::indexOf(graph.residues.nodes.nodeAdjacencies[residueA], residueB);
            const std::vector<size_t>& edges = graph.residues.nodes.edgeAdjacencies[residueA];
            if (adjacency < edges.size())
            {
                size_t edgeIndex = edges[adjacency];
                const std::array<std::vector<bool>, 2>& ignoredAtoms = residueAtomsCloseToEdge[edgeIndex];
                bool order = !(graph.residues.edges.nodeAdjacencies[edgeIndex][0] == residueA);
                nonIgnoredA = util::vectorNot(ignoredAtoms[order]);
                nonIgnoredB = util::vectorNot(ignoredAtoms[!order]);
            }
            for (size_t n : util::boolsToValues(atomsA, nonIgnoredA))
            {
                for (size_t k : util::boolsToValues(atomsB, nonIgnoredB))
                {
                    Element elementA = data.atoms.elements[n];
                    Element elementB = data.atoms.elements[k];
                    if (elementRadii.hasValue[elementA] && elementRadii.hasValue[elementB])
                    {
                        Sphere a = bounds.atoms[n];
                        Sphere b = bounds.atoms[k];
                        double cutoff = std::max(0.0, a.radius + b.radius - overlapTolerance);
                        double sqDist = squaredDistance(a.center, b.center);
                        if (sqDist < cutoff * cutoff)
                        {
                            PotentialFactor factor = potentialFactor(potential, elementA, elementB);
                            double lennardJones = lennardJonesPotential(factor, sqDist);
                            result.push_back({n, k, (a.radius + b.radius) - std::sqrt(sqDist), lennardJones});
                        }
                    }
                }
            }
        }

        return result;
    };
    std::function<std::vector<std::string>(const AtomContact&)> contactToRow = [&](const AtomContact& contact)
    {
        size_t atomA = contact.atomA;
        size_t atomB = contact.atomB;
        size_t residueA = data.indices.atomResidue[atomA];
        size_t residueB = data.indices.atomResidue[atomB];
        return std::vector<std::string> {
            data.residues.chainIds[residueA],
            data.residues.names[residueA],
            std::to_string(data.residues.numbers[residueA]),
            data.atoms.names[atomA],
            data.residues.chainIds[residueB],
            data.residues.names[residueB],
            std::to_string(data.residues.numbers[residueB]),
            data.atoms.names[atomB],
            std::to_string(contact.overlap),
            std::to_string(contact.potential)};
    };
    std::vector<std::string> contactHeader {"", "", "", "", "", "", "", "", "overlap", "lennard-jones potential"};

    util::TextTable textTable {{}, {}};

    switch (mode)
    {
        case MOLECULES:
            {
                textTable = {
                    moleculeHeader, util::vectorMap(moleculeRow, util::indexVector(data.indices.moleculeCount))};
                break;
            }
        case RESIDUES:
            {
                textTable = {residueHeader, util::vectorMap(residueRow, util::indexVector(data.indices.residueCount))};
                break;
            }
        case ATOMS:
            {
                textTable = {
                    atomHeader,
                    util::vectorMap(atomRow, util::indicesOfElement(data.atoms.elements, Element::Unknown))};
                break;
            }
        case CONTACTS:
            {
                size_t residueCount = data.indices.residueCount;
                std::vector<AtomContact> contacts;
                for (size_t n = 0; n < residueCount; n++)
                {
                    for (size_t k = n + 1; k < residueCount; k++)
                    {
                        util::insertInto(contacts, contactRows(n, k));
                    }
                }
                std::function<bool(const AtomContact&, const AtomContact&)> highestOverlap =
                    [](const AtomContact& a, const AtomContact& b) { return a.overlap > b.overlap; };
                std::function<bool(const AtomContact&, const AtomContact&)> highestPotential =
                    [](const AtomContact& a, const AtomContact& b) { return a.potential > b.potential; };
                if (order == OVERLAP)
                {
                    contacts = util::sortedBy(highestOverlap, contacts);
                }
                else if (order == POTENTIAL)
                {
                    contacts = util::sortedBy(highestPotential, contacts);
                }
                textTable = {contactHeader, util::vectorMap(contactToRow, contacts)};
                break;
            }
    }
    switch (format)
    {
        case CSV:
            util::toCsv(std::cout, ",", textTable);
            break;
        case TXT:
            util::toTxt(std::cout, {textTable});
            break;
        case HTML:
            util::toHtml(std::cout, {textTable});
            break;
    }
    return 0;
}
