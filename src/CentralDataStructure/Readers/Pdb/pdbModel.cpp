#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"

#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CentralDataStructure/Editors/amberMdPrep.hpp" //all preprocessing should move to here.
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbChain.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Writers/print.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

void pdb::readAssembly(PdbData& data, size_t assemblyId, cds::Assembly& assembly, std::stringstream& stream_block)
{
    int currentModelNumber = 1;
    std::string previousResidueId = "InitialValue";
    std::string line;
    while (getline(stream_block, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if (recordName == "MODEL")
        {
            try
            {
                currentModelNumber = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(10, 4)));
            }
            catch (...)
            {
                gmml::log(
                    __LINE__,
                    __FILE__,
                    gmml::WAR,
                    "Model issue: this ain't an int: " + codeUtils::RemoveWhiteSpace(line.substr(10, 4)));
                currentModelNumber = 1; // Seems like a reasonable default.
            }
            data.assemblies.numbers[assemblyId] = currentModelNumber;
            assembly.setNumber(currentModelNumber);
        }
        else if ((recordName == "ATOM") || (recordName == "HETATM"))
        { // Gimme everything with the same chain, can be everything with no chain.
            // Function that will read from stringstream until chain ID changes or TER or just not ATOM/HETATM
            std::string chainId = extractChainId(line);
            std::stringstream singleChainSection = extractSingleChainFromRecordSection(stream_block, line, chainId);
            size_t moleculeId = data.indices.moleculeCount;
            data.indices.moleculeAssembly.push_back(assemblyId);
            data.indices.moleculeCount++;
            data.molecules.chainIds.push_back(chainId);
            data.molecules.residueOrder.push_back({});
            std::unique_ptr<cds::Molecule> molecule = std::make_unique<cds::Molecule>();
            cds::Molecule* mol = assembly.addMolecule(std::move(molecule));
            data.objects.molecules.push_back(mol);
            readChain(data, moleculeId, singleChainSection);
        }
        else if (recordName == "ENDMDL")
        { // Only happens when reading "modelsAsCoordinates", now read the rest of the entry as extra coords for
          // MODEL 1.
            gmml::log(__LINE__, __FILE__, gmml::INF, "PdbFile being read in as trajectory");
            extractCoordinatesFromModel(data, stream_block, line);
        }
    }
    return;
}

std::string pdb::extractChainId(const std::string& line)
{ // serialNumber can overrun into position 12 in input.
    int shift = checkShiftFromSerialNumberOverrun(line);
    return codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
}

std::stringstream pdb::extractSingleChainFromRecordSection(
    std::stringstream& stream_block, std::string line, const std::string& initialChainID)
{
    std::streampos previousLinePosition = stream_block.tellg(); // Save current line position
    std::stringstream singleChainSection;
    std::string chainID = initialChainID;
    std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
    while ((chainID == initialChainID) && (recordName != "TER"))
    {
        singleChainSection << line << std::endl;
        previousLinePosition = stream_block.tellg(); // Save current line position.
        if (!std::getline(stream_block, line))
        {
            break; // // If we hit the end, time to leave.
        }
        chainID = extractChainId(line);
        recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
    }
    stream_block.seekg(
        previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    return singleChainSection;
}

void pdb::extractCoordinatesFromModel(PdbData& data, std::stringstream& stream_block, std::string line)
{
    const int iPdbLineLength = 80; // repeat for now, fix later
    gmml::log(__LINE__, __FILE__, gmml::INF, "Section to extract coordinates from is\n" + stream_block.str());
    size_t atomCount = data.indices.atomCount;
    if (atomCount == 0)
    {
        std::string message = "No atoms available when extracting coords from multiple models";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw message;
    }
    const std::vector<cds::Coordinate>& initial = data.atoms.coordinates;
    data.trajectory.coordinates = {initial};
    std::vector<cds::Coordinate> current = initial;
    size_t atomId = 0;
    while ((std::getline(stream_block, line)))
    {
        expandLine(line, iPdbLineLength);
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if (recordName == "ATOM" || recordName == "HETATM")
        {
            if (atomId >= atomCount)
            {
                std::string message = "Trajectory model contains more atoms than initial model";
                gmml::log(__LINE__, __FILE__, gmml::ERR, message);
                throw message;
            }
            current[atomId] = checkShiftsAndExtractCoordinate(line);
            atomId++;
        }
        if (recordName == "ENDMDL")
        { // reset to read next set of coords
            atomId = 0;
            data.trajectory.coordinates.push_back(current);
            current = initial;
        }
    }
    return;
}

void pdb::Write(const PdbData& data, const std::vector<std::vector<size_t>>& moleculeResidues, std::ostream& stream)
{
    std::function<bool(const size_t&)> atomAlive = [&](size_t n) { return data.indices.atomAlive[n]; };
    std::function<std::string(const size_t&)> elementString = [&](size_t n)
    { return data.atoms.names[n].empty() ? "" : data.atoms.names[n].substr(0, 1); };
    std::function<std::string(const size_t&)> truncatedResidueNames = [&](size_t n)
    { return codeUtils::truncate(3, data.residues.names[n]); };

    assembly::Graph graph = cds::createAssemblyGraph(data.indices, data.atomGraph);
    cds::PdbFileResidueData residueData {
        data.residues.numbers,
        codeUtils::vectorMap(truncatedResidueNames, codeUtils::indexVector(data.indices.residueCount)),
        data.residues.chainIds,
        data.residues.insertionCodes};
    cds::PdbFileAtomData atomData {
        data.atoms.coordinates,
        data.atoms.numbers,
        data.atoms.names,
        codeUtils::vectorMap(elementString, codeUtils::indexVector(data.indices.atomCount)),
        data.atoms.recordNames,
        data.atoms.occupancies,
        data.atoms.temperatureFactors};
    cds::PdbFileFormat format;
    cds::PdbFileData writerData {format, {}, residueData, atomData};
    for (auto& residueIds : moleculeResidues)
    {
        for (size_t residueId : residueIds)
        {
            cds::writeAssemblyToPdb(
                stream, graph, {{residueId}}, {{data.residues.hasTerCard[residueId]}}, {}, writerData);
        }
        if (!residueIds.empty())
        { // Sometimes you get empty chains after things have been deleted I guess.
            stream << "TER\n";
        }
    }
}
