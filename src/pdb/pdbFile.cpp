#include "include/pdb/pdbFile.hpp"

#include "include/CentralDataStructure/assembly.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/pdb/pdbFileWriter.hpp"
#include "include/pdb/pdbFunctions.hpp"
#include "include/pdb/pdbModel.hpp"
#include "include/pdb/pdbRecords.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/logging.hpp"
#include "include/util/parsing.hpp"
#include "include/util/strings.hpp"

#include <fstream>
#include <iomanip> // setprecision setw
#include <ostream>

namespace gmml
{
    namespace pdb
    {
        namespace
        {
            void readConectRow(PdbData& data, const std::string& line)
            {
                auto indexOfAtom = [&](const std::string& str)
                {
                    std::optional<int> parsed = util::parseInt(str);
                    if (!parsed.has_value())
                    {
                        throw std::runtime_error("Error: could not parse conect id: " + str);
                    }
                    uint number = parsed.value();
                    size_t index = util::indexOf(data.atoms.numbers, number);
                    if (index == data.atoms.names.size())
                    {
                        throw std::runtime_error("Error: conect row atom id not found: " + str);
                    }
                    return index;
                };
                std::vector<std::string> split = util::split(line, ' ');
                size_t first = indexOfAtom(split[1]);
                for (size_t n = 2; n < split.size(); n++)
                {
                    addBond(data, first, indexOfAtom(split[n]));
                }
            }

            // Goes through a section of the PDB file that contains the same header section. E.g. HEADER.
            // If the header changes, it goes back to the previous line. I wanted the out while loop to trigger the new
            // line. This means I don't have to check recordName between if statement and can have else if.
            std::stringstream extractHomogenousRecordSection(
                std::istream& pdbFileStream, std::string& line, std::string recordName)
            {
                std::stringstream recordSection;
                expandLine(line, iPdbLineLength);
                recordSection << line << std::endl;
                std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
                std::string previousName = recordName;
                while (std::getline(pdbFileStream, line))
                {
                    expandLine(line, iPdbLineLength);
                    recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                    if (recordName != "ANISOU") // Do nothing for ANISOU
                    {
                        if (recordName == previousName)
                        {
                            recordSection << line << std::endl;
                            previousName = recordName;
                            previousLinePosition = pdbFileStream.tellg(); // Save current line position.
                        }
                        else
                        {
                            break;
                        }
                    } // Do nothing for ANISOU
                }
                pdbFileStream.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER
                                                           // and found TITLE.
                return recordSection;
            }

            // Initializers used by constructors
            // Should extract all lines that start with the strings in recordNames.
            // Returns when it hits a line that does not start with one of those records.
            std::stringstream extractHeterogenousRecordSection(
                std::istream& pdbFileStream, std::string& line, const std::vector<std::string> recordNames)
            {
                std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
                std::stringstream recordSection;
                std::string recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                while (util::contains(recordNames, recordName))
                {
                    if (recordName != "ANISOU") // Do nothing for ANISOU
                    {
                        std::stringstream partialRecordSection =
                            extractHomogenousRecordSection(pdbFileStream, line, recordName);
                        recordSection << partialRecordSection.str();
                        previousLinePosition = pdbFileStream.tellg(); // Save current line position.
                    }
                    if (!std::getline(pdbFileStream, line)) // If we hit the end
                    {
                        break; // Time to leave.
                    }
                    recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                }
                pdbFileStream.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER
                                                           // and found TITLE.
                return recordSection;
            }

            void parseInFileStream(PdbFile& file, std::istream& pdbFileStream, const ReaderOptions& options)
            {
                PdbData& data = file.data;
                size_t assemblyId = 0;
                for (std::string line; std::getline(pdbFileStream, line);)
                {
                    expandLine(line, iPdbLineLength);
                    std::string recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                    std::vector<std::string> coordSectionCards {"MODEL", "ATOM", "ANISOU", "TER", "HETATM"};
                    if (options.inputType == modelsAsCoordinates)
                    { // want to pass in the whole block to assembly so it can read the extra coords
                        coordSectionCards.push_back("ENDMDL");
                    }
                    std::vector<std::string> databaseCards {"DBREF", "DBREF1", "DBREF2"};
                    if (util::contains(coordSectionCards, recordName))
                    {
                        std::stringstream recordSection =
                            extractHeterogenousRecordSection(pdbFileStream, line, coordSectionCards);
                        Assembly& assembly = file.assemblies.emplace_back(Assembly());
                        data.indices.assemblyCount++;
                        data.assemblies.numbers.push_back(assemblyId + 1);
                        readAssembly(data, assemblyId, assembly, recordSection);
                        assemblyId++;
                    }
                    else if (recordName == "HEADER")
                    {
                        std::stringstream recordSection =
                            extractHomogenousRecordSection(pdbFileStream, line, recordName);
                        file.headerRecord = readHeaderRecord(recordSection);
                    }
                    else if (recordName == "TITLE")
                    {
                        std::stringstream recordSection =
                            extractHomogenousRecordSection(pdbFileStream, line, recordName);
                        file.titleRecord = readTitleRecord(recordSection);
                    }
                    else if (recordName == "AUTHOR")
                    {
                        std::stringstream recordSection =
                            extractHomogenousRecordSection(pdbFileStream, line, recordName);
                        file.authorRecord = readAuthorRecord(recordSection);
                    }
                    else if (recordName == "JRNL")
                    {
                        std::stringstream recordSection =
                            extractHomogenousRecordSection(pdbFileStream, line, recordName);
                        file.journalRecord = readJournalRecord(recordSection);
                    }
                    else if (recordName == "REMARK")
                    {
                        std::stringstream recordSection =
                            extractHomogenousRecordSection(pdbFileStream, line, recordName);
                        file.remarkRecord = readRemarkRecord(recordSection);
                    }
                    else if (util::contains(databaseCards, recordName))
                    {
                        std::stringstream databaseSection =
                            extractHeterogenousRecordSection(pdbFileStream, line, databaseCards);
                        while (getline(databaseSection, line))
                        {
                            file.databaseReferences.push_back(readDatabaseReference(line));
                        }
                    }
                    else if (recordName == "CONECT")
                    {
                        if (options.readConectRows)
                        {
                            readConectRow(data, line);
                        }
                        else
                        {
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::WAR,
                                "Reading pdbfile that contains CONECT records. We ignore these due to the potential "
                                "for "
                                "overruns.");
                        }
                    }
                }
                util::log(__LINE__, __FILE__, util::INF, "PdbFile Constructor Complete Captain");
                data.objects.assemblies = getAssemblies(file);
            }
        } // namespace

        PdbFile toPdbFile(const std::string& pdbFilePath, const ReaderOptions& options)
        {
            PdbFile result;
            result.inFilePath = pdbFilePath;
            std::ifstream pdbFileStream(pdbFilePath);
            if (pdbFileStream.fail())
            {
                util::log(__LINE__, __FILE__, util::ERR, "Could not open this file: " + pdbFilePath);
                throw std::runtime_error("PdbFile constructor could not open this file: " + pdbFilePath);
            }
            util::log(__LINE__, __FILE__, util::INF, "File opened: " + pdbFilePath + ". Ready to parse!");
            parseInFileStream(result, pdbFileStream, options);
            util::log(__LINE__, __FILE__, util::INF, "Finished parsing " + pdbFilePath);
            return result;
        }

        std::vector<Assembly*> getAssemblies(PdbFile& file)
        {
            std::vector<Assembly*> result;
            result.reserve(file.assemblies.size());
            for (auto& a : file.assemblies)
            {
                result.push_back(&a);
            }
            return result;
        }

        void write(PdbFile& file, const std::string outName)
        {
            util::writeToFile(outName, [&](std::ostream& stream) { write(file, stream); });
        }

        void write(PdbFile& file, std::ostream& out)
        {
            PdbData& data = file.data;
            write(file.headerRecord, out);
            write(file.titleRecord, out);
            write(file.authorRecord, out);
            write(file.journalRecord, out);
            write(file.remarkRecord, out);
            for (auto& dbref : file.databaseReferences)
            {
                write(dbref, out);
            }
            for (size_t n = 0; n < data.indices.assemblyCount; n++)
            {
                if (data.indices.assemblyCount > 1)
                {
                    out << "MODEL " << std::right << std::setw(4) << data.assemblies.numbers[n] << "\n";
                }
                std::vector<size_t> moleculeIds = assemblyMolecules(data.indices, n);
                pdb::write(data, util::indicesToValues(data.molecules.residueOrder, moleculeIds), out);
                if (data.indices.assemblyCount > 1)
                {
                    out << "ENDMDL\n";
                }
            }
            theEnd(out);
        }
    } // namespace pdb
} // namespace gmml
