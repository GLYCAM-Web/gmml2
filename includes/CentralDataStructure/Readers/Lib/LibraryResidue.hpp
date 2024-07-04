#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYRESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYRESIDUE_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CentralDataStructure/Readers/Lib/LibraryAtom.hpp"
#include <sstream>

namespace lib
{
    class LibraryResidue : public cds::Residue
    {
      public:
        LibraryResidue(std::stringstream& residueStream, const std::string& name);

      private:
        //  I don't think we ever use the box
        //    double box_angle_;
        //    double box_length_;
        //    double box_width_;
        //    double box_height_;
        int head_atom_index_ = constants::iNotSet;
        int tail_atom_index_ = constants::iNotSet;
    };
} // namespace lib
#endif
