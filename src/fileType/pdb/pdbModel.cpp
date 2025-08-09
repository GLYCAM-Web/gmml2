#include "include/fileType/pdb/pdbModel.hpp"

#include "include/CentralDataStructure/assembly.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/fileType/pdb/bondByDistance.hpp"
#include "include/fileType/pdb/pdbChain.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/fileType/pdb/pdbFileWriter.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

namespace gmml
{
    namespace pdb
    {
        void readAssembly(PdbData& data, size_t assemblyId, Assembly& assembly, std::stringstream& stream_block)
        {
            int currentModelNumber = 1;
            std::string previousResidueId = "InitialValue";
            std::string line;
            while (getline(stream_block, line))
            {
                std::string recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                if (recordName == "MODEL")
                {
                    try
                    {
                        currentModelNumber = std::stoi(util::RemoveWhiteSpace(line.substr(10, 4)));
                    }
                    catch (...)
                    {
                        util::log(
                            __LINE__,
                            __FILE__,
                            util::WAR,
                            "Model issue: this ain't an int: " + util::RemoveWhiteSpace(line.substr(10, 4)));
                        currentModelNumber = 1; // Seems like a reasonable default.
                    }
                    data.assemblies.numbers[assemblyId] = currentModelNumber;
                    assembly.setNumber(currentModelNumber);
                }
                else if ((recordName == "ATOM") || (recordName == "HETATM"))
                { // Gimme everything with the same chain, can be everything with no chain.
                    // Function that will read from stringstream until chain ID changes or TER or just not ATOM/HETATM
                    std::string chainId = extractChainId(line);
                    std::stringstream singleChainSection =
                        extractSingleChainFromRecordSection(stream_block, line, chainId);
                    size_t moleculeId = data.indices.moleculeCount;
                    data.indices.moleculeAssembly.push_back(assemblyId);
                    data.indices.moleculeCount++;
                    data.molecules.chainIds.push_back(chainId);
                    data.molecules.residueOrder.push_back({});
                    std::unique_ptr<Molecule> molecule = std::make_unique<Molecule>();
                    Molecule* mol = assembly.addMolecule(std::move(molecule));
                    data.objects.molecules.push_back(mol);
                    readChain(data, moleculeId, singleChainSection);
                }
                else if (recordName == "ENDMDL")
                { // Only happens when reading "modelsAsCoordinates", now read the rest of the entry as extra coords for
                  // MODEL 1.
                    util::log(__LINE__, __FILE__, util::INF, "PdbFile being read in as trajectory");
                    extractCoordinatesFromModel(data, stream_block, line);
                }
            }
            return;
        }

        std::string extractChainId(const std::string& line)
        { // serialNumber can overrun into position 12 in input.
            int shift = checkShiftFromSerialNumberOverrun(line);
            return util::RemoveWhiteSpace(line.substr(21 + shift, 1));
        }

        std::stringstream extractSingleChainFromRecordSection(
            std::stringstream& stream_block, std::string line, const std::string& initialChainID)
        {
            std::streampos previousLinePosition = stream_block.tellg(); // Save current line position
            std::stringstream singleChainSection;
            std::string chainID = initialChainID;
            std::string recordName = util::RemoveWhiteSpace(line.substr(0, 6));
            while ((chainID == initialChainID) && (recordName != "TER"))
            {
                singleChainSection << line << std::endl;
                previousLinePosition = stream_block.tellg(); // Save current line position.
                if (!std::getline(stream_block, line))
                {
                    break; // // If we hit the end, time to leave.
                }
                chainID = extractChainId(line);
                recordName = util::RemoveWhiteSpace(line.substr(0, 6));
            }
            stream_block.seekg(
                previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
            return singleChainSection;
        }

        void extractCoordinatesFromModel(PdbData& data, std::stringstream& stream_block, std::string line)
        {
            const int iPdbLineLength = 80; // repeat for now, fix later
            util::log(__LINE__, __FILE__, util::INF, "Section to extract coordinates from is\n" + stream_block.str());
            size_t atomCount = data.indices.atomCount;
            if (atomCount == 0)
            {
                std::string message = "No atoms available when extracting coords from multiple models";
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw message;
            }
            const std::vector<Coordinate>& initial = data.atoms.coordinates;
            data.trajectory.coordinates = {initial};
            std::vector<Coordinate> current = initial;
            size_t atomId = 0;
            while ((std::getline(stream_block, line)))
            {
                expandLine(line, iPdbLineLength);
                std::string recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                if (recordName == "ATOM" || recordName == "HETATM")
                {
                    if (atomId >= atomCount)
                    {
                        std::string message = "Trajectory model contains more atoms than initial model";
                        util::log(__LINE__, __FILE__, util::ERR, message);
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

        void write(const PdbData& data, const std::vector<std::vector<size_t>>& moleculeResidues, std::ostream& stream)
        {
            std::function<bool(const size_t&)> atomAlive = [&](size_t n) { return data.indices.atomAlive[n]; };
            std::function<std::string(const size_t&)> elementString = [&](size_t n)
            { return data.atoms.names[n].empty() ? "" : data.atoms.names[n].substr(0, 1); };
            std::function<std::string(const size_t&)> truncatedResidueNames = [&](size_t n)
            { return util::truncate(3, data.residues.names[n]); };

            assembly::Graph graph = createAssemblyGraph(data.indices, data.atomGraph);
            PdbFileResidueData residueData {
                data.residues.numbers,
                util::vectorMap(truncatedResidueNames, util::indexVector(data.indices.residueCount)),
                data.residues.chainIds,
                data.residues.insertionCodes};
            PdbFileAtomData atomData {
                data.atoms.coordinates,
                data.atoms.numbers,
                data.atoms.names,
                util::vectorMap(elementString, util::indexVector(data.indices.atomCount)),
                data.atoms.recordNames,
                data.atoms.occupancies,
                data.atoms.temperatureFactors};
            PdbFileFormat format;
            PdbFileData writerData {format, {}, residueData, atomData};
            for (auto& residueIds : moleculeResidues)
            {
                for (size_t residueId : residueIds)
                {
                    writeAssemblyToPdb(
                        stream, graph, {{residueId}}, {{data.residues.hasTerCard[residueId]}}, {}, writerData);
                }
                if (!residueIds.empty())
                { // Sometimes you get empty chains after things have been deleted I guess.
                    stream << "TER\n";
                }
            }
        }
    } // namespace pdb
} // namespace gmml
