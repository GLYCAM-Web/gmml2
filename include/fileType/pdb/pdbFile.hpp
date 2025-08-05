#ifndef INCLUDE_FILETYPE_PDB_PDBFILE_HPP
#define INCLUDE_FILETYPE_PDB_PDBFILE_HPP
// ToDo split preprocessor into separate function.
// ToDo Warning about gaps being 1 residue.
// ToDo make more direct queries here instead of giving out HeaderRecord etc.
// ToDo ACE/NME between residues with same number but an insertion code.
#include "include/CentralDataStructure/assembly.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbRecordTypes.hpp"

#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
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

        struct PdbFile
        {
            std::string inFilePath = "";
            HeaderRecord headerRecord;
            TitleRecord titleRecord;
            AuthorRecord authorRecord;
            JournalRecord journalRecord;
            RemarkRecord remarkRecord;
            std::vector<DatabaseReference> databaseReferences;
            std::vector<Assembly> assemblies;
            PdbData data;
        };

        PdbFile toPdbFile(const std::string& pdbFilePath, const ReaderOptions& options);

        inline PdbFile toPdbFile(const std::string& pdbFilePath, const InputType pdbFileType)
        {
            return toPdbFile(pdbFilePath, {pdbFileType, false});
        }

        std::vector<Assembly*> getAssemblies(PdbFile& file);
        void write(PdbFile& file, const std::string outName);
        void write(PdbFile& file, std::ostream& stream);
    } // namespace pdb
} // namespace gmml

#endif
