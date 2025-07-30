#ifndef INCLUDE_PDB_PDBRECORDTYPES_HPP
#define INCLUDE_PDB_PDBRECORDTYPES_HPP

#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        struct AuthorRecord
        {
            std::string recordName; /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::string author;     /*!< Author that appears in KEYWORD record of a pdb file >*/
        };

        struct DatabaseReference
        {
            std::string record_name_;  /*!< Record name of database reference which is the first column of each line of
                                        the card >*/
            std::string id_code_;      /*!< ID code of this entry >*/
            std::string chain_id_;     /*!< Chain ID of the database reference >*/
            int seq_begin_;            /*!< Initial sequence number of the PDB sequence segment, right justified >*/
            std::string insert_begin_; /*!< Initial insertion code of the PDB sequence segment >*/
            int seq_end_;              /*!< Ending sequence number of the PDB sequence segment, right justified >*/
            std::string insert_end_;   /*!< Ending insertion code of the PDB sequence segment >*/
            std::string database_;     /*!< Sequence database name >*/
            std::string db_accession_; /*!< Sequence database accession code >*/
            std::string db_id_code_;   /*!< Sequence  database identification code >*/
            int db_seq_begin_;         /*!< Initial sequence number of the database segment >*/
            std::string db_ins_beg_; /*!< Insertion code of initial residue of the segment, if PDB is the reference >*/
            int db_seq_end_;         /*!< Ending sequence number of the database segment >*/
            std::string
                db_ins_end_; /*!< Insertion code of the ending residue of the segment, if PDB is the reference >*/
        };

        struct HeaderRecord
        {
            std::string recordName;     /*!< Record name of headr card in a pdb file >*/
            std::string classification; /*!< Classification of the pdb file >*/
            std::string depositionDate; /*!< Date of deposition >*/
            std::string identifierCode; /*!< Identifier code of the pdb file >*/
        };

        struct JournalRecord
        {
            std::string recordName =
                ""; /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::vector<std::string> authors; /*!< List of authors that appear in Journal record of a pdb file >*/
            std::string title;                /*!< Title that appears in Journal record of a pdb file >*/
            std::vector<std::string> editors; /*!< List of editors that appear in Journal record of a pdb file >*/
            std::string reference;            /*!< Reference that appears in Journal record of a pdb file >*/
            std::string publisher;            /*!< Publisher that appears in Journal record of a pdb file >*/
            std::vector<std::string>
                referenceNums; /*!< List of reference numbers that appear in Journal record of a pdb file >*/
            std::string pmid;  /*!< Pub Med ID number that appears in Journal record of a pdb file >*/
            std::string doi;   /*!< DOI number that appears in Journal record of a pdb file >*/
            std::string text;  /*!< Text in a Journal Section >*/
        };

        struct RemarkRecord
        {
            double resolution; /*!< Resolution of PDB >*/
            double b_factor;   /*!< B Factor of PDB >*/
        };

        struct TitleRecord
        {
            std::string name;  /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::string title; /*!< Title that appears in TITLE record of a pdb file >*/
        };
    } // namespace pdb
} // namespace gmml

#endif
