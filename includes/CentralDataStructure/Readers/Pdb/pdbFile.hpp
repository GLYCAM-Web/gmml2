#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFILE_HPP
// ToDo split preprocessor into separate function.
// ToDo Warning about gaps being 1 residue.
// ToDo make more direct queries here instead of giving out HeaderRecord etc.
// ToDo ACE/NME between residues with same number but an insertion code.
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/authorRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/databaseReferenceRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/headerRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/journalRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/remarkRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/titleRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/assembly.hpp"

#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace pdb
{
    const int iPdbLineLength = 80;

    enum InputType
    {
        modelsAsMolecules,
        modelsAsCoordinates,
    };

    struct ReaderOptions
    {
        InputType inputType;
        bool readConectRows;
    };

    class PdbFile
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        PdbFile();
        PdbFile(const std::string& pdbFilePath, const InputType pdbFileType = modelsAsMolecules);
        PdbFile(const std::string& pdbFilePath, const ReaderOptions& options);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline std::string GetInputFilePath() const { return inFilePath_; }

        // ToDo These should be private and whatever info they give out should be directly queryable here.
        inline const HeaderRecord& GetHeaderRecord() const { return headerRecord_; }

        inline const TitleRecord& GetTitleRecord() const { return titleRecord_; }

        inline const AuthorRecord& GetAuthorRecord() const { return authorRecord_; }

        inline const JournalRecord& GetJournalRecord() const { return journalRecord_; }

        inline const RemarkRecord& GetRemarkRecord() const { return remarkRecord_; }

        inline std::vector<cds::Assembly*> getAssemblies()
        {
            std::vector<cds::Assembly*> result;
            result.reserve(assemblies_.size());
            for (auto& a : assemblies_)
            {
                result.push_back(&a);
            }
            return result;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string GetUniprotIDs() const;
        const float& GetResolution() const;
        const float& GetBFactor() const;
        PreprocessorInformation PreProcess(
            const cdsParameters::ParameterManager& parameterManager, PreprocessorOptions options);
        //////////////////////////////////////////////////////////
        //                        DISPLAY                       //
        //////////////////////////////////////////////////////////
        void Write(const std::string outName);
        void Write(std::ostream& stream);

      private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline const std::vector<DatabaseReference>& GetDatabaseReferences() const { return databaseReferences_; }

        //////////////////////////////////////////////////////////
        //                        ATTRIBUTES                    //
        //////////////////////////////////////////////////////////
      public:
        std::string inFilePath_ = "";
        HeaderRecord headerRecord_; // SWIG wants the
        TitleRecord titleRecord_;
        AuthorRecord authorRecord_;
        JournalRecord journalRecord_;
        RemarkRecord remarkRecord_;
        std::vector<DatabaseReference> databaseReferences_;
        std::vector<cds::Assembly> assemblies_;
        PdbData data;
    };
} // namespace pdb
#endif
