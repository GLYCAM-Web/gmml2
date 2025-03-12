#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_SECTIONCLASSES_TITLERECORD_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_SECTIONCLASSES_TITLERECORD_HPP

#include <string>
#include <sstream>
#include <ostream>

namespace pdb
{
    class TitleRecord
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        TitleRecord(std::string name = "TITLE", std::string title = "DEFAULT");
        TitleRecord(std::stringstream& stream_block);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        const std::string& GetRecordName() const;
        const std::string& GetTitle() const;
        //////////////////////////////////////////////////////////
        //                        FUNCTIONS                     //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        void Print(std::ostream& out) const;
        void Write(std::ostream& stream) const;

      private:
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void SetName(const std::string name);
        void SetTitle(const std::string title);
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        std::string name_  = ""; /*!< Record name which appears in the first column of each line in a pdb file >*/
        std::string title_ = ""; /*!< Title that appears in TITLE record of a pdb file >*/
    };
} // namespace pdb

#endif
