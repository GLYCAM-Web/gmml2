#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyIndices.hpp"
#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/structuredFiles.hpp"
#include "includes/version.h"

#include <cmath>
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

    using codeUtils::ArgReq;
    using codeUtils::ArgType;
    // keeping the same order between the enum and strings lets us assign format by index in the vector
    std::vector<std::string> knownFormats {"csv", "txt", "html"};
    std::vector<std::string> modes {"molecules", "residues", "atoms", "contacts"};
    std::vector<std::string> orders {"index", "overlap", "potential"};
    std::vector<codeUtils::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed, FILENAME, "", ' ', "filename"},
        {ArgReq::optional, ArgType::flag, HELP, "help", 'h', ""},
        {ArgReq::optional, ArgType::flag, VERSION, "version", 'v', ""},
        {ArgReq::optional, ArgType::option, FORMAT, "format", 'f', codeUtils::join("|", knownFormats)},
        {ArgReq::optional, ArgType::option, MODE, "mode", 'm', codeUtils::join("|", modes)},
        {ArgReq::optional, ArgType::option, ORDER, "order", 'o', codeUtils::join("|", orders)},
        {ArgReq::optional, ArgType::option, OVERLAPTOLERANCE, "tolerance", 't', "value"}
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
            std::cout << "PDB Analyzer & GMML version " << GMML_VERSION << "\n";
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
    OUTPUT_FORMAT format      = TXT;
    OUTPUT_MODE mode          = CONTACTS;
    OUTPUT_ORDER order        = INDEX;
    double overlapTolerance   = 0.0;
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
            case ARGUMENTS::MODE:
                {
                    size_t index = codeUtils::indexOf(modes, arg.value);
                    if (index >= modes.size())
                    {
                        throw std::runtime_error("Unknown mode: '" + arg.value + "', expected one of " +
                                                 codeUtils::join(", ", modes));
                    }
                    mode = OUTPUT_MODE(index);
                    break;
                }
            case ARGUMENTS::ORDER:
                {
                    size_t index = codeUtils::indexOf(orders, arg.value);
                    if (index >= orders.size())
                    {
                        throw std::runtime_error("Unknown order: '" + arg.value + "', expected one of " +
                                                 codeUtils::join(", ", orders));
                    }
                    order = OUTPUT_ORDER(index);
                    break;
                }
            case ARGUMENTS::OVERLAPTOLERANCE:
                {
                    std::optional<double> opt = codeUtils::parseDouble(arg.value);
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
    using MolecularMetadata::Element;
    pdb::PdbFile inputFile(inputFileName);
    pdb::PdbData& data = inputFile.data;
    pdb::bondAtomsAndResiduesByDistance(data);
    assembly::Graph graph                           = cds::createAssemblyGraph(data.indices, data.atomGraph);
    auto residueAtomsCloseToEdge                    = assembly::atomsCloseToResidueEdges(graph);
    const codeUtils::SparseVector<double> atomRadii = MolecularMetadata::vanDerWaalsRadii();
    const MolecularMetadata::PotentialTable potential =
        MolecularMetadata::potentialTable(atomRadii, MolecularMetadata::foundElements(data.atoms.elements));
    std::function<cds::Sphere(const size_t&)> toAtomBounds = [&](size_t atomId)
    {
        return cds::Sphere {atomRadii.values[data.atoms.elements[atomId]], data.atoms.coordinates[atomId]};
    };

    std::vector<cds::Sphere> atomBounds =
        codeUtils::vectorMap(toAtomBounds, codeUtils::indexVector(data.indices.atomCount));
    const assembly::Bounds bounds = cds::toAssemblyBounds(graph, atomBounds);

    std::vector<std::string> residueTypeNames(cds::ResidueTypeCount, "");
    residueTypeNames[cds::ResidueType::Aglycone]            = "aglycone";
    residueTypeNames[cds::ResidueType::Deoxy]               = "deoxy";
    residueTypeNames[cds::ResidueType::Derivative]          = "derivative";
    residueTypeNames[cds::ResidueType::Protein]             = "protein";
    residueTypeNames[cds::ResidueType::ProteinCappingGroup] = "protein cap";
    residueTypeNames[cds::ResidueType::Solvent]             = "solvent";
    residueTypeNames[cds::ResidueType::Sugar]               = "sugar";
    residueTypeNames[cds::ResidueType::Undefined]           = "unknown";

    std::function<std::vector<std::string>(const size_t&)> moleculeRow = [&](size_t moleculeId)
    {
        std::string chain = data.molecules.chainIds[moleculeId];
        std::vector<bool> includedTypes(cds::ResidueTypeCount, false);
        for (cds::ResidueType type :
             codeUtils::indicesToValues(data.residues.types, moleculeResidues(data.indices, moleculeId)))
        {
            includedTypes[type] = true;
        }
        std::string types        = codeUtils::join(", ", codeUtils::boolsToValues(residueTypeNames, includedTypes));
        std::string residueCount = std::to_string(assembly::moleculeResidues(data.indices, moleculeId).size());
        std::string atomCount    = std::to_string(assembly::moleculeAtoms(data.indices, moleculeId).size());
        return std::vector<std::string> {chain, residueCount, atomCount, types};
    };
    std::vector<std::string> moleculeHeader {"chain", "residues", "atoms", "residue types"};

    std::function<std::vector<std::string>(const size_t&)> residueRow = [&](size_t residueId)
    {
        std::string chain     = data.residues.chainIds[residueId];
        std::string name      = data.residues.names[residueId];
        std::string number    = std::to_string(data.residues.numbers[residueId]);
        std::string type      = residueTypeNames[data.residues.types[residueId]];
        std::string atomCount = std::to_string(assembly::residueAtoms(data.indices, residueId).size());
        return std::vector<std::string> {chain, name, number, type, atomCount};
    };
    std::vector<std::string> residueHeader {"chain", "name", "number", "type", "atom count"};

    std::function<std::vector<std::string>(const size_t&)> atomRow = [&](size_t atomId)
    {
        size_t residueId          = data.indices.atomResidue[atomId];
        std::string chain         = data.residues.chainIds[residueId];
        std::string residueName   = data.residues.names[residueId];
        std::string residueNumber = std::to_string(data.residues.numbers[residueId]);
        std::string name          = data.atoms.names[atomId];
        std::string number        = std::to_string(data.atoms.numbers[atomId]);
        std::string element       = MolecularMetadata::elementName(data.atoms.elements[atomId]);
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
        bool bothProtein = (data.residues.types[residueA] == cds::ResidueType::Protein) &&
                           (data.residues.types[residueB] == cds::ResidueType::Protein);
        bool eitherWater =
            codeUtils::contains({data.residues.names[residueA], data.residues.names[residueB]}, std::string("HOH"));
        if (!bothProtein && !eitherWater &&
            cds::spheresOverlap(overlapTolerance, bounds.residues[residueA], bounds.residues[residueB]))
        {
            const std::vector<size_t>& atomsA = graph.residues.nodes.constituents[residueA];
            const std::vector<size_t>& atomsB = graph.residues.nodes.constituents[residueB];
            std::vector<bool> nonIgnoredA(atomsA.size(), true);
            std::vector<bool> nonIgnoredB(atomsB.size(), true);
            size_t adjacency = codeUtils::indexOf(graph.residues.nodes.nodeAdjacencies[residueA], residueB);
            const std::vector<size_t>& edges = graph.residues.nodes.edgeAdjacencies[residueA];
            if (adjacency < edges.size())
            {
                size_t edgeIndex                                     = edges[adjacency];
                const std::array<std::vector<bool>, 2>& ignoredAtoms = residueAtomsCloseToEdge[edgeIndex];
                bool order  = !(graph.residues.edges.nodeAdjacencies[edgeIndex][0] == residueA);
                nonIgnoredA = codeUtils::vectorNot(ignoredAtoms[order]);
                nonIgnoredB = codeUtils::vectorNot(ignoredAtoms[!order]);
            }
            for (size_t n : codeUtils::boolsToValues(atomsA, nonIgnoredA))
            {
                for (size_t k : codeUtils::boolsToValues(atomsB, nonIgnoredB))
                {
                    Element elementA = data.atoms.elements[n];
                    Element elementB = data.atoms.elements[k];
                    if (atomRadii.hasValue[elementA] && atomRadii.hasValue[elementB])
                    {
                        cds::Sphere a = bounds.atoms[n];
                        cds::Sphere b = bounds.atoms[k];
                        double cutoff = std::max(0.0, a.radius + b.radius - overlapTolerance);
                        double sqDist = cds::squaredDistance(a.center, b.center);
                        if (sqDist < cutoff * cutoff)
                        {
                            MolecularMetadata::PotentialFactor factor =
                                MolecularMetadata::potentialFactor(potential, elementA, elementB);
                            double lennardJones = MolecularMetadata::lennardJonesPotential(factor, sqDist);
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
        size_t atomA    = contact.atomA;
        size_t atomB    = contact.atomB;
        size_t residueA = data.indices.atomResidue[atomA];
        size_t residueB = data.indices.atomResidue[atomB];
        return std::vector<std::string> {data.residues.chainIds[residueA],
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

    codeUtils::TextTable textTable {{}, {}};

    switch (mode)
    {
        case MOLECULES:
            {
                textTable = {moleculeHeader,
                             codeUtils::vectorMap(moleculeRow, codeUtils::indexVector(data.indices.moleculeCount))};
                break;
            }
        case RESIDUES:
            {
                textTable = {residueHeader,
                             codeUtils::vectorMap(residueRow, codeUtils::indexVector(data.indices.residueCount))};
                break;
            }
        case ATOMS:
            {
                textTable = {atomHeader, codeUtils::vectorMap(
                                             atomRow, codeUtils::indicesOfElement(
                                                          data.atoms.elements, MolecularMetadata::Element::Unknown))};
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
                        codeUtils::insertInto(contacts, contactRows(n, k));
                    }
                }
                std::function<bool(const AtomContact&, const AtomContact&)> highestOverlap =
                    [](const AtomContact& a, const AtomContact& b)
                {
                    return a.overlap > b.overlap;
                };
                std::function<bool(const AtomContact&, const AtomContact&)> highestPotential =
                    [](const AtomContact& a, const AtomContact& b)
                {
                    return a.potential > b.potential;
                };
                if (order == OVERLAP)
                {
                    contacts = codeUtils::sortedBy(highestOverlap, contacts);
                }
                else if (order == POTENTIAL)
                {
                    contacts = codeUtils::sortedBy(highestPotential, contacts);
                }
                textTable = {contactHeader, codeUtils::vectorMap(contactToRow, contacts)};
                break;
            }
    }
    switch (format)
    {
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
