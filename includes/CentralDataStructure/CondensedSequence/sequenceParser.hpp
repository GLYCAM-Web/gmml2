#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"

#include <memory>
#include <string>
#include <vector>

namespace cdsCondensedSequence
{
    void parseSequence(std::vector<std::unique_ptr<ParsedResidue>>& residues, std::string sequence);
} // namespace cdsCondensedSequence
#endif
