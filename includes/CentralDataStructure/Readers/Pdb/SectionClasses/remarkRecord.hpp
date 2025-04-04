#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_SECTIONCLASSES_REMARKRECORD_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_SECTIONCLASSES_REMARKRECORD_HPP

#include <string>
#include <sstream>
#include <ostream>

namespace pdb
{
    class RemarkRecord
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        RemarkRecord();
        RemarkRecord(std::stringstream& stream_block);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline const float& GetResolution() const
        {
            return resolution_;
        }

        inline const float& GetBFactor() const
        {
            return b_factor_;
        }

        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        void Print(std::ostream& out) const;
        void Write(std::ostream& stream) const;

      private:
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void SetResolution(const float resolution);
        void SetBFactor(const float b_factor);
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        float resolution_ = 0.0; /*!< Resolution of PDB >*/
        float b_factor_   = 0.0; /*!< B Factor of PDB >*/
    };
} // namespace pdb

#endif
