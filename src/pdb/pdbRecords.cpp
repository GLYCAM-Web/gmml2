#include "include/pdb/pdbRecords.hpp"

#include "include/pdb/pdbRecordTypes.hpp"
#include "include/util/constants.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        AuthorRecord readAuthorRecord(std::stringstream& stream_block)
        {
            std::string line;
            bool is_record_name_set = false;
            std::stringstream ss;
            getline(stream_block, line);
            std::string temp = line;
            std::string recordName;
            while (!util::Trim(temp).empty())
            {
                if (!is_record_name_set)
                {
                    recordName = line.substr(0, 6);
                    util::Trim(recordName);
                    is_record_name_set = true;
                }
                ss << line.substr(10, 70);

                getline(stream_block, line);
                temp = line;
            }
            return {recordName, ss.str()};
        }

        void write(const AuthorRecord& record, std::ostream& stream)
        {
            const int MAX_AUTHOR_LENGTH_IN_LINE = 70;
            const std::string& author = record.author;
            const std::string& recordName = record.recordName;
            stream << std::left << std::setw(6) << recordName << std::left << std::setw(2) << " ";
            if ((int)author.length() > MAX_AUTHOR_LENGTH_IN_LINE)
            {
                stream << std::left << std::setw(70) << author.substr(0, MAX_AUTHOR_LENGTH_IN_LINE) << std::endl;

                int counter = ceil((double)(author.length()) / MAX_AUTHOR_LENGTH_IN_LINE);
                for (int i = 2; i <= counter; i++)
                {
                    if (i != counter)
                    {
                        stream << std::left << std::setw(6) << recordName << std::left << std::setw(2) << " "
                               << std::right << std::setw(2) << i << std::left << std::setw(70)
                               << author.substr(MAX_AUTHOR_LENGTH_IN_LINE * (i - 1), MAX_AUTHOR_LENGTH_IN_LINE)
                               << std::endl;
                    }
                    else
                    {
                        stream << std::left << std::setw(6) << recordName << std::left << std::setw(2) << " "
                               << std::right << std::setw(2) << i << std::left << std::setw(70)
                               << author.substr(
                                      MAX_AUTHOR_LENGTH_IN_LINE * (i - 1),
                                      author.length() - MAX_AUTHOR_LENGTH_IN_LINE * (i - 1))
                               << std::endl;
                    }
                }
            }
            else
            {
                stream << std::right << std::setw(2) << " " << std::left << std::setw(70) << author << std::endl;
            }
        }

        DatabaseReference readDatabaseReference(const std::string& line)
        {
            DatabaseReference result;
            std::string temp = line;
            if (!util::Trim(temp).empty())
            {
                result.record_name_ = line.substr(0, 6);
                result.id_code_ = line.substr(7, 4);
                result.chain_id_ = line.substr(12, 1);
                result.seq_begin_ = util::from_string<int>(line.substr(14, 4));
                result.insert_begin_ = line.substr(18, 1);
                result.seq_end_ = util::from_string<int>(line.substr(20, 4));
                result.insert_end_ = line.substr(24, 1);
                result.database_ = line.substr(26, 6);

                if (result.record_name_ == "DBREF ")
                {
                    result.db_accession_ = line.substr(33, 8);
                    result.db_id_code_ = line.substr(42, 12);
                    result.db_seq_begin_ = util::from_string<int>(line.substr(55, 5));
                    result.db_ins_beg_ = line.substr(60, 1);
                    result.db_seq_end_ = util::from_string<int>(line.substr(62, 5));
                    result.db_ins_end_ = line.substr(67, 1);
                }
                else if (result.record_name_ == "DBREF1")
                {
                    std::size_t position = line.find("DBREF2");
                    if (position != std::string::npos)
                    {
                        std::string dbref2 = line.substr(position, 68);
                        result.db_accession_ = dbref2.substr(18, 22);
                        result.db_id_code_ = line.substr(47, 15);
                        result.db_seq_begin_ = util::from_string<int>(dbref2.substr(45, 10));
                        result.db_ins_beg_ = ' ';
                        result.db_seq_end_ = util::from_string<int>(line.substr(57, 10));
                        result.db_ins_end_ = ' ';
                    }
                }
            }
            return result;
        }

        void write(const DatabaseReference& record, std::ostream& stream)
        {
            const std::string& recordName = record.record_name_;
            const std::string& idCode = record.id_code_;
            const std::string& dbIdCode = record.db_id_code_;
            const std::string& chainId = record.chain_id_;
            int seqBegin = record.seq_begin_;
            int seqEnd = record.seq_end_;
            int dbSeqBegin = record.db_seq_begin_;
            int dbSeqEnd = record.db_seq_end_;
            const std::string& insertBegin = record.insert_begin_;
            const std::string& insertEnd = record.insert_end_;
            const std::string& dbInsertBegin = record.db_ins_beg_;
            const std::string& dbInsertEnd = record.db_ins_end_;
            const std::string& database = record.database_;
            const std::string& dbAccess = record.db_accession_;
            if (recordName == "DBREF ")
            {
                stream << std::left << std::setw(6) << recordName << std::left << std::setw(1) << " " << std::left
                       << std::setw(4) << idCode << std::left << std::setw(1) << " " << std::left << std::setw(1)
                       << chainId << std::left << std::setw(1) << " " << std::right << std::setw(4) << seqBegin
                       << std::right << std::setw(1) << insertBegin << std::left << std::setw(1) << " " << std::right
                       << std::setw(4) << seqEnd << std::right << std::setw(1) << insertEnd << std::left << std::setw(1)
                       << " " << std::left << std::setw(6) << database << std::left << std::setw(1) << " " << std::left
                       << std::setw(8) << dbAccess << std::left << std::setw(1) << " " << std::left << std::setw(12)
                       << dbIdCode << std::left << std::setw(1) << " " << std::right << std::setw(5) << dbSeqBegin
                       << std::right << std::setw(1) << dbInsertBegin << std::left << std::setw(1) << " " << std::right
                       << std::setw(5) << dbSeqEnd << std::right << std::setw(1) << dbInsertEnd << std::endl;
            }
            else if (recordName == "DBREF1")
            {
                stream << std::left << std::setw(6) << recordName << std::left << std::setw(1) << " " << std::left
                       << std::setw(4) << idCode << std::left << std::setw(1) << " " << std::left << std::setw(1)
                       << chainId << std::left << std::setw(1) << " " << std::right << std::setw(4) << seqBegin
                       << std::right << std::setw(1) << insertBegin << std::left << std::setw(1) << " " << std::right
                       << std::setw(4) << seqEnd << std::right << std::setw(1) << insertEnd << std::left << std::setw(1)
                       << " " << std::left << std::setw(6) << database << std::left << std::setw(16) << " " << std::left
                       << std::setw(15) << dbIdCode << std::endl
                       << std::left << std::setw(6) << "DBREF2" << std::left << std::setw(1) << " " << std::left
                       << std::setw(4) << idCode << std::left << std::setw(1) << " " << std::left << std::setw(1)
                       << chainId << std::left << std::setw(6) << " " << std::left << std::setw(22) << dbAccess
                       << std::left << std::setw(5) << " " << std::right << std::setw(10) << dbSeqBegin << std::left
                       << std::setw(2) << " " << std::right << std::setw(10) << dbSeqEnd << std::endl;
            }
        }

        HeaderRecord readHeaderRecord(std::stringstream& stream_block)
        {
            std::string recordName;
            std::string classification;
            std::string depositionDate;
            std::string identifierCode;
            std::string line;
            getline(stream_block, line);
            std::string temp = line;
            while (!util::Trim(temp).empty())
            {
                recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                classification = util::RemoveWhiteSpace(line.substr(10, 40));
                depositionDate = util::RemoveWhiteSpace(line.substr(50, 9));
                identifierCode = util::RemoveWhiteSpace(line.substr(62, 4));
                getline(stream_block, line);
                temp = line;
            }
            return {recordName, classification, depositionDate, identifierCode};
        }

        void write(const HeaderRecord& record, std::ostream& stream)
        {
            stream << std::left << std::setw(6) << record.recordName << std::left << std::setw(4) << " " << std::left
                   << std::setw(40) << record.classification << std::left << std::setw(9) << record.depositionDate
                   << std::left << std::setw(3) << " " << std::right << std::setw(4) << record.identifierCode
                   << std::left << std::setw(14) << " " << std::endl;
        }

        JournalRecord readJournalRecord(std::stringstream& stream_block)
        {
            std::string line;
            bool is_record_name_set = false;
            bool is_title_started = false;
            bool is_reference_started = false;
            bool is_publisher_started = false;
            std::string recordName;
            std::string text;
            std::string title;
            std::string reference;
            std::string publisher;
            std::string pmid;
            std::string doi;
            std::vector<std::string> authors;
            std::vector<std::string> editors;
            std::vector<std::string> referenceNums;
            getline(stream_block, line);
            std::string temp = line;
            while (!util::Trim(temp).empty())
            {
                if (!is_record_name_set)
                {
                    recordName = line.substr(0, 6);
                    util::Trim(recordName);
                    is_record_name_set = true;
                }
                text.append(line.substr(12, 67));
                std::string subrecord = line.substr(12, 4);
                util::Trim(subrecord);
                if (subrecord == "AUTH")
                {
                    std::size_t start_position = 19;
                    std::size_t end_position = line.find(",");
                    if (end_position != std::string::npos)
                    {
                        std::string new_author = line.substr(start_position, end_position - start_position);
                        authors.push_back(new_author);
                        start_position = end_position;
                        end_position = line.find(",", start_position + 1);
                    }
                    else
                    {
                        std::string new_author = line.substr(start_position, 79 - start_position);
                        authors.push_back(new_author);
                    }
                }
                else if (subrecord == "TITL")
                {
                    if (!is_title_started)
                    {
                        title = line.substr(19, 59);
                        util::Trim(title);
                        is_title_started = true;
                    }
                    else
                    {
                        std::string s = line.substr(19, 59);
                        util::Trim(s);
                        title.append(" ");
                        title.append(s);
                        util::Trim(title);
                    }
                }
                else if (subrecord == "EDIT")
                {
                    std::size_t start_position = 19;
                    std::size_t end_position = line.find(",");
                    while (end_position != std::string::npos)
                    {
                        std::string new_editor = line.substr(start_position, end_position - start_position);
                        editors.push_back(new_editor);
                        start_position = end_position;
                        end_position = line.find(",", start_position + 1);
                    }
                }
                else if (subrecord == "REF")
                {
                    if (!is_reference_started)
                    {
                        reference = line.substr(19, 59);
                        util::Trim(reference);
                        is_reference_started = true;
                    }
                    else
                    {
                        std::string s = line.substr(19, 59);
                        util::Trim(s);
                        reference.append(" ");
                        reference.append(s);
                        util::Trim(reference);
                    }
                }
                else if (subrecord == "PUBL")
                {
                    if (!is_publisher_started)
                    {
                        publisher = line.substr(19, 59);
                        util::Trim(publisher);
                        is_publisher_started = true;
                    }
                    else
                    {
                        std::string publ = line.substr(19, 59);
                        util::Trim(publ);
                        publisher.append(" ");
                        publisher.append(publ);
                    }
                }
                else if (subrecord == "REFN")
                {
                    std::string refn = line.substr(35, 29);
                    util::Trim(refn);
                    referenceNums.push_back(refn);
                }
                else if (subrecord == "PMID")
                {
                    std::string new_pmid = line.substr(19, 59);
                    util::Trim(new_pmid);
                    pmid.append(new_pmid);
                }
                else if (subrecord == "DOI")
                {
                    std::string new_doi = line.substr(19, 59);
                    doi += new_doi;
                    util::Trim(doi);
                    util::removeMultipleSpaces(doi);
                    // I know this is weird to Trim spaces then add one space back, but it was/is necessary. DT
                    doi += " ";
                }
                getline(stream_block, line);
                temp = line;
            }
            return {recordName, authors, title, editors, reference, publisher, referenceNums, pmid, doi, text};
        }

        void write(const JournalRecord& record, std::ostream& stream)
        {
            const int MAX_JOURNAL_LENGTH_IN_LINE = 67;
            const std::string& text = record.text;
            const std::string& recordName = record.recordName;
            stream << std::left << std::setw(6) << recordName << std::left << std::setw(6) << " ";
            if ((int)text.length() > MAX_JOURNAL_LENGTH_IN_LINE)
            {
                stream << std::left << std::setw(67) << text.substr(0, MAX_JOURNAL_LENGTH_IN_LINE) << std::endl;

                int counter = ceil((double)(text.length()) / MAX_JOURNAL_LENGTH_IN_LINE);
                for (int i = 2; i <= counter; i++)
                {
                    if (i != counter)
                    {
                        stream << std::left << std::setw(6) << recordName << std::left << std::setw(6) << " "
                               << std::left << std::setw(67)
                               << text.substr(MAX_JOURNAL_LENGTH_IN_LINE * (i - 1), MAX_JOURNAL_LENGTH_IN_LINE)
                               << std::endl;
                    }
                    else
                    {
                        stream << std::left << std::setw(6) << recordName << std::left << std::setw(6) << " "
                               << std::left << std::setw(67)
                               << text.substr(
                                      MAX_JOURNAL_LENGTH_IN_LINE * (i - 1),
                                      text.length() - MAX_JOURNAL_LENGTH_IN_LINE * (i - 1))
                               << std::endl;
                    }
                }
            }
            else
            {
                stream << std::right << std::setw(6) << " " << std::left << std::setw(67) << text << std::endl;
            }
        }

        RemarkRecord readRemarkRecord(std::stringstream& stream_block)
        {
            double resolution = constants::dNotSet;
            double bfactor = constants::dNotSet;
            std::string line;
            getline(stream_block, line);
            std::string temp = line;
            while (!util::Trim(temp).empty())
            {
                if (line.find("REMARK") != std::string::npos)
                {
                    if (line.find("2 RESOLUTION.") != std::string::npos)
                    {
                        std::string tmp_resolution = line.substr(23, 7);
                        util::Trim(tmp_resolution);
                        try
                        {
                            resolution = std::stof(tmp_resolution);
                        }
                        catch (const std::invalid_argument& error)
                        {
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::ERR,
                                "RESOLUTION is not a valid float value. Value:\t" + tmp_resolution);
                        }
                    }
                    if (line.find("MEAN B VALUE") != std::string::npos)
                    {
                        int start = line.find(":") + 1;
                        std::string tmp_b_factor = line.substr(start, 80 - start);
                        util::Trim(tmp_b_factor);
                        try
                        {
                            bfactor = std::stof(tmp_b_factor);
                        }
                        catch (const std::invalid_argument& error)
                        {
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::ERR,
                                "MEAN B VALUE is not a valid float value. Value:\t" + tmp_b_factor);
                        }
                    }
                }
                getline(stream_block, line);
                temp = line;
            }
            return {resolution, bfactor};
        }

        void write(const RemarkRecord& record, std::ostream& stream)
        {
            stream << "REMARK   2\n";
            stream << "REMARK   2 RESOLUTION.  " << record.resolution << " ANGSTROMS.\n";
            stream << "REMARK   3\n";
            stream << "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : " << record.b_factor << "\n";
        }

        TitleRecord readTitleRecord(std::stringstream& stream_block)
        {
            std::string line;
            bool is_name_set = false;
            std::stringstream ss;
            getline(stream_block, line);
            std::string temp = line;
            std::string name;
            while (!util::Trim(temp).empty())
            {
                if (!is_name_set)
                {
                    name = line.substr(0, 6);
                    util::Trim(name);
                    is_name_set = true;
                }
                ss << line.substr(10, 70);

                getline(stream_block, line);
                temp = line;
            }
            std::string title = ss.str();
            util::Trim(title);
            util::removeMultipleSpaces(title);
            return {name, title};
        }

        void write(const TitleRecord& record, std::ostream& stream)
        { // OG: I just copied this...
            const int MAX_TITLE_LENGTH_IN_LINE = 70;
            const std::string& name = record.name;
            const std::string& title = record.title;
            stream << std::left << std::setw(6) << name << std::left << std::setw(2) << " ";
            if ((int)title.length() > MAX_TITLE_LENGTH_IN_LINE)
            {
                stream << std::right << std::setw(2) << " " << std::left << std::setw(70)
                       << title.substr(0, MAX_TITLE_LENGTH_IN_LINE) << std::endl;

                int counter = ceil((double)(title.length()) / MAX_TITLE_LENGTH_IN_LINE);
                for (int i = 2; i <= counter; i++)
                {
                    if (i != counter)
                    {
                        stream << std::left << std::setw(6) << name << std::left << std::setw(2) << " " << std::right
                               << std::setw(2) << i << std::left << std::setw(70)
                               << title.substr(MAX_TITLE_LENGTH_IN_LINE * (i - 1), MAX_TITLE_LENGTH_IN_LINE)
                               << std::endl;
                    }
                    else
                    {
                        stream << std::left << std::setw(6) << name << std::left << std::setw(2) << " " << std::right
                               << std::setw(2) << i << std::left << std::setw(70)
                               << title.substr(
                                      MAX_TITLE_LENGTH_IN_LINE * (i - 1),
                                      title.length() - MAX_TITLE_LENGTH_IN_LINE * (i - 1))
                               << std::endl;
                    }
                }
            }
            else
            {
                stream << std::right << std::setw(2) << " " << std::left << std::setw(70) << title << std::endl;
            }
        }
    } // namespace pdb
} // namespace gmml
