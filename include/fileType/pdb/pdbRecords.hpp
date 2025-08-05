#ifndef INCLUDE_FILETYPE_PDB_PDBRECORDS_HPP
#define INCLUDE_FILETYPE_PDB_PDBRECORDS_HPP

#include "include/fileType/pdb/pdbRecordTypes.hpp"

#include <ostream>
#include <sstream>
#include <string>

namespace gmml
{
    namespace pdb
    {
        AuthorRecord readAuthorRecord(std::stringstream& stream_block);
        DatabaseReference readDatabaseReference(const std::string& line);
        HeaderRecord readHeaderRecord(std::stringstream& stream_block);
        JournalRecord readJournalRecord(std::stringstream& stream_block);
        RemarkRecord readRemarkRecord(std::stringstream& stream_block);
        TitleRecord readTitleRecord(std::stringstream& stream_block);

        void write(const AuthorRecord& record, std::ostream& stream);
        void write(const DatabaseReference& record, std::ostream& stream);
        void write(const HeaderRecord& record, std::ostream& stream);
        void write(const JournalRecord& record, std::ostream& stream);
        void write(const RemarkRecord& record, std::ostream& stream);
        void write(const TitleRecord& record, std::ostream& stream);
    } // namespace pdb
} // namespace gmml

#endif
