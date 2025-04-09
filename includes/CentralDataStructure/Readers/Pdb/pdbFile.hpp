#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFILE_HPP
// ToDo split preprocessor into separate function.
// ToDo Warning about gaps being 1 residue.
// ToDo make more direct queries here instead of giving out HeaderRecord etc.
// ToDo ACE/NME between residues with same number but an insertion code.
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/authorRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/databaseReferenceRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/headerRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/journalRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/remarkRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/titleRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include <string>
#include <istream>
#include <ostream>
#include <vector>

namespace pdb
{
    const int iPdbLineLength = 80;

    enum InputType
    {
        modelsAsMolecules,
        modelsAsCoordinates,
    };

    class PdbFile
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        PdbFile();
        PdbFile(const std::string& pdbFilePath, const InputType pdbFileType = modelsAsMolecules);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline std::string GetInputFilePath() const
        {
            return inFilePath_;
        }

        // ToDo These should be private and whatever info they give out should be directly queryable here.
        inline const HeaderRecord& GetHeaderRecord() const
        {
            return headerRecord_;
        }

        inline const TitleRecord& GetTitleRecord() const
        {
            return titleRecord_;
        }

        inline const AuthorRecord& GetAuthorRecord() const
        {
            return authorRecord_;
        }

        inline const JournalRecord& GetJournalRecord() const
        {
            return journalRecord_;
        }

        inline const RemarkRecord& GetRemarkRecord() const
        {
            return remarkRecord_;
        }

        inline const std::vector<PdbModel>& getAssemblies() const
        {
            return assemblies_;
        }

        inline std::vector<PdbModel>& mutableAssemblies()
        {
            return assemblies_;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string GetUniprotIDs() const;
        const float& GetResolution() const;
        const float& GetBFactor() const;
        PreprocessorInformation PreProcess(const cdsParameters::ParameterManager& parameterManager,
                                           PreprocessorOptions options);
        //////////////////////////////////////////////////////////
        //                        DISPLAY                       //
        //////////////////////////////////////////////////////////
        void Write(const std::string outName) const;
        void Write(std::ostream& stream) const;

      private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline const std::vector<DatabaseReference>& GetDatabaseReferences() const
        {
            return databaseReferences_;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void ParseInFileStream(std::istream& pdbFileStream, const InputType pdbFileType = modelsAsMolecules);
        std::stringstream ExtractHeterogenousRecordSection(std::istream& pdbFileStream, std::string& line,
                                                           const std::vector<std::string> recordNames);
        std::stringstream ExtractHomogenousRecordSection(std::istream& pdbFileStream, std::string& line,
                                                         std::string previousName);
        //////////////////////////////////////////////////////////
        //                        ATTRIBUTES                    //
        //////////////////////////////////////////////////////////
        std::string inFilePath_ = "";
        HeaderRecord headerRecord_; // SWIG wants the
        TitleRecord titleRecord_;
        AuthorRecord authorRecord_;
        JournalRecord journalRecord_;
        RemarkRecord remarkRecord_;
        std::vector<DatabaseReference> databaseReferences_;
        std::vector<PdbModel> assemblies_;
    };
} // namespace pdb
#endif
