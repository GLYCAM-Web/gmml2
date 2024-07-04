#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_SECTIONCLASSES_CONECTRECORD_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_SECTIONCLASSES_CONECTRECORD_HPP

#include <string>
#include <iostream>

#include "includes/CentralDataStructure/atom.hpp"

namespace pdb
{
    class PdbModel;

    class ConectRecord
    {
        // This is passed in as serial numbers, but they can/will change so store pointers to the atoms that are bonded
        // instead. When asked to print or output to file, get the serial numbers via the pointers.
      public:
        //////////////////////////////////////////////////////////
        //                    CONSTRUCTOR                       //
        //////////////////////////////////////////////////////////
        ConectRecord(std::string& line, PdbModel& pdbModel);
        ConectRecord(std::vector<const cds::Atom*> atomRecords);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        void Print(std::ostream& out = std::cerr) const;
        void Write(std::ostream& stream) const;

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        std::vector<const cds::Atom*> atomRecordPtrs_;
    };
} // namespace pdb

#endif
