#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulation.hpp"

#include <string>

namespace cdsCondensedSequence
{
    AbstractSequence parseSequence(std::string sequence);

    inline SequenceData parseAndReorder(const std::string& sequence)
    {
        return reordered(instantiate(parseSequence(sequence)));
    }

} // namespace cdsCondensedSequence
#endif
