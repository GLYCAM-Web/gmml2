#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreProcess.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <fstream>
#include <ostream>
#include <iomanip> // setprecision setw

using pdb::PdbFile;

namespace
{
    using pdb::PdbData;
    using pdb::ReaderOptions;

    void readConectRow(PdbData& data, const std::string& line)
    {
        auto indexOfAtom = [&](const std::string& str)
        {
            std::optional<int> parsed = codeUtils::parseInt(str);
            if (!parsed.has_value())
            {
                throw std::runtime_error("Error: could not parse conect id: " + str);
            }
            uint number  = parsed.value();
            size_t index = codeUtils::indexOf(data.atoms.numbers, number);
            if (index == data.atoms.names.size())
            {
                throw std::runtime_error("Error: conect row atom id not found: " + str);
            }
            return index;
        };
        std::vector<std::string> split = codeUtils::split(line, ' ');
        size_t first                   = indexOfAtom(split[1]);
        for (size_t n = 2; n < split.size(); n++)
        {
            addBond(data, first, indexOfAtom(split[n]));
        }
    }

    // Goes through a section of the PDB file that contains the same header section. E.g. HEADER.
    // If the header changes, it goes back to the previous line. I wanted the out while loop to trigger the new line.
    // This means I don't have to check recordName between if statement and can have else if.
    std::stringstream extractHomogenousRecordSection(std::istream& pdbFileStream, std::string& line,
                                                     std::string recordName)
    {
        std::stringstream recordSection;
        pdb::expandLine(line, pdb::iPdbLineLength);
        recordSection << line << std::endl;
        std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
        std::string previousName            = recordName;
        while ((std::getline(pdbFileStream, line)))
        {
            pdb::expandLine(line, pdb::iPdbLineLength);
            recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
            if (recordName != "ANISOU") // Do nothing for ANISOU
            {
                if (recordName == previousName)
                {
                    recordSection << line << std::endl;
                    previousName         = recordName;
                    previousLinePosition = pdbFileStream.tellg(); // Save current line position.
                }
                else
                {
                    break;
                }
            } // Do nothing for ANISOU
        }
        pdbFileStream.seekg(
            previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
        // std::cout << "At end. Returning this record section:\n" << recordSection.str() << "\nEND RECORD SECTION\n";
        return recordSection;
    }

    // Initializers used by constructors
    // Should extract all lines that start with the strings in recordNames.
    // Returns when it hits a line that does not start with one of those records.
    std::stringstream extractHeterogenousRecordSection(std::istream& pdbFileStream, std::string& line,
                                                       const std::vector<std::string> recordNames)
    {
        std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
        std::stringstream recordSection;
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        while (codeUtils::contains(recordNames, recordName))
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
            recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        }
        pdbFileStream.seekg(
            previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
        //    std::cout << "Returning this hetero section:\n" << heteroRecordSection.str() << "\nThe End of Hetero
        //    Section.\n";
        return recordSection;
    }

    void parseInFileStream(PdbFile& file, std::istream& pdbFileStream, const ReaderOptions& options)
    {
        PdbData& data     = file.data;
        size_t assemblyId = 0;
        for (std::string line; std::getline(pdbFileStream, line);)
        {
            pdb::expandLine(line, pdb::iPdbLineLength);
            std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
            std::vector<std::string> coordSectionCards {"MODEL", "ATOM", "ANISOU", "TER", "HETATM"};
            if (options.inputType == pdb::modelsAsCoordinates)
            { // want to pass in the whole block to assembly so it can read the extra coords
                coordSectionCards.push_back("ENDMDL");
            }
            std::vector<std::string> databaseCards {"DBREF", "DBREF1", "DBREF2"};
            if (codeUtils::contains(coordSectionCards, recordName))
            {
                std::stringstream recordSection =
                    extractHeterogenousRecordSection(pdbFileStream, line, coordSectionCards);
                cds::Assembly& assembly = file.assemblies_.emplace_back(cds::Assembly());
                data.indices.assemblyCount++;
                data.assemblies.numbers.push_back(assemblyId + 1);
                readAssembly(data, assemblyId, assembly, recordSection);
                assemblyId++;
            }
            else if (recordName == "HEADER")
            {
                std::stringstream recordSection = extractHomogenousRecordSection(pdbFileStream, line, recordName);
                file.headerRecord_              = pdb::HeaderRecord(recordSection);
            }
            else if (recordName == "TITLE")
            {
                std::stringstream recordSection = extractHomogenousRecordSection(pdbFileStream, line, recordName);
                file.titleRecord_               = pdb::TitleRecord(recordSection);
            }
            else if (recordName == "AUTHOR")
            {
                std::stringstream recordSection = extractHomogenousRecordSection(pdbFileStream, line, recordName);
                file.authorRecord_              = pdb::AuthorRecord(recordSection);
            }
            else if (recordName == "JRNL")
            {
                std::stringstream recordSection = extractHomogenousRecordSection(pdbFileStream, line, recordName);
                file.journalRecord_             = pdb::JournalRecord(recordSection);
            }
            else if (recordName == "REMARK")
            {
                std::stringstream recordSection = extractHomogenousRecordSection(pdbFileStream, line, recordName);
                file.remarkRecord_              = pdb::RemarkRecord(recordSection);
            }
            else if (codeUtils::contains(databaseCards, recordName))
            {
                std::stringstream databaseSection =
                    extractHeterogenousRecordSection(pdbFileStream, line, databaseCards);
                while (getline(databaseSection, line))
                {
                    file.databaseReferences_.emplace_back(line);
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
                    gmml::log(__LINE__, __FILE__, gmml::WAR,
                              "Reading pdbfile that contains CONECT records. We ignore these due to the potential for "
                              "overruns.");
                }
            }
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "PdbFile Constructor Complete Captain");
        data.objects.assemblies = file.getAssemblies();
    }
} // namespace

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFile::PdbFile()
{
    inFilePath_ = "GMML-Generated";
}

PdbFile::PdbFile(const std::string& pdbFilePath, const InputType pdbFileType)
    : PdbFile(pdbFilePath, {pdbFileType, false})
{}

PdbFile::PdbFile(const std::string& pdbFilePath, const ReaderOptions& options) : inFilePath_(pdbFilePath)
{
    std::ifstream pdbFileStream(pdbFilePath);
    if (pdbFileStream.fail())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Could not open this file: " + pdbFilePath);
        throw std::runtime_error("PdbFile constructor could not open this file: " + pdbFilePath);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "File opened: " + pdbFilePath + ". Ready to parse!");
    parseInFileStream(*this, pdbFileStream, options);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished parsing " + pdbFilePath);
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
std::string PdbFile::GetUniprotIDs() const
{
    std::string UniprotIDs = "";
    for (auto& databaseReference : this->GetDatabaseReferences())
    {
        UniprotIDs += databaseReference.GetUniprotID();
    }
    return UniprotIDs;
}

const float& PdbFile::GetResolution() const
{
    return this->GetRemarkRecord().GetResolution();
}

const float& PdbFile::GetBFactor() const
{
    return this->GetRemarkRecord().GetBFactor();
}

pdb::PreprocessorInformation PdbFile::PreProcess(const cdsParameters::ParameterManager& parameterManager,
                                                 PreprocessorOptions inputOptions)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Preprocesssing has begun");
    pdb::PreprocessorInformation ppInfo;
    for (size_t assemblyId = 0; assemblyId < assemblies_.size();
         assemblyId++) // Now we do all, but maybe user can select at some point.
    {
        cds::Assembly* assembly = &assemblies_[assemblyId];
        preProcessCysResidues(data, assemblyId, ppInfo);
        preProcessHisResidues(data, assemblyId, ppInfo, inputOptions);
        preProcessChainTerminals(data, assemblyId, ppInfo, inputOptions);
        preProcessGapsUsingDistance(data, assemblyId, ppInfo, inputOptions);
        preProcessMissingUnrecognized(data, assemblyId, ppInfo, parameterManager);
        cdsParameters::setAtomChargesForResidues(parameterManager, assembly->getResidues());
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Preprocessing completed");
    return ppInfo;
}

void PdbFile::Write(const std::string outName)
{
    codeUtils::writeToFile(outName,
                           [&](std::ostream& stream)
                           {
                               this->Write(stream);
                           });
}

void PdbFile::Write(std::ostream& out)
{
    this->GetHeaderRecord().Write(out);
    this->GetTitleRecord().Write(out);
    this->GetAuthorRecord().Write(out);
    this->GetJournalRecord().Write(out);
    this->GetRemarkRecord().Write(out);
    for (auto& dbref : this->GetDatabaseReferences())
    {
        dbref.Write(out);
    }
    for (size_t n = 0; n < data.indices.assemblyCount; n++)
    {
        if (data.indices.assemblyCount > 1)
        {
            out << "MODEL " << std::right << std::setw(4) << data.assemblies.numbers[n] << "\n";
        }
        std::vector<size_t> moleculeIds = assemblyMolecules(data.indices, n);
        pdb::Write(data, codeUtils::indicesToValues(data.molecules.residueOrder, moleculeIds), out);
        if (data.indices.assemblyCount > 1)
        {
            out << "ENDMDL\n";
        }
    }
    cds::theEnd(out);
}
