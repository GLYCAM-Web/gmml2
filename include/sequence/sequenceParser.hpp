#ifndef INCLUDES_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP
#define INCLUDES_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP

#include "include/sequence/sequenceManipulation.hpp"
#include "include/sequence/sequenceTypes.hpp"

#include <string>

namespace gmml
{
    namespace sequence
    {
        AbstractSequence parseSequence(std::string sequence);

        inline SequenceData parseAndReorder(const std::string& sequence)
        {
            return reordered(instantiate(parseSequence(sequence)));
        }
    } // namespace sequence
} // namespace gmml

#endif
