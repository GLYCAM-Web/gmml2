#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP

#include "includes/CentralDataStructure/molecule.hpp"

#include <string>

namespace cdsCondensedSequence
{
    void parseSequence(cds::Molecule* molecule, std::string sequence);
} // namespace cdsCondensedSequence
#endif
