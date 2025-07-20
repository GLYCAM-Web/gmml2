#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"

#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulation.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequencePrinter.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>

CondensedSequence::Sequence::Sequence(std::string condensedSequence)
{
    cdsCondensedSequence::AbstractSequence abstract = cdsCondensedSequence::parseSequence(condensedSequence);
    cdsCondensedSequence::SequenceData sequence = cdsCondensedSequence::instantiate(abstract);
    this->setInterpretedSequence(
        cdsCondensedSequence::printSequence(cdsCondensedSequence::defaultConfigUnsorted(), sequence));
    cdsCondensedSequence::SequenceData reordered = cdsCondensedSequence::reordered(sequence);
    this->setInputSequence(condensedSequence);
    this->setIndexOrdered(cdsCondensedSequence::printSequence(cdsCondensedSequence::defaultConfig(), reordered));
    this->setIndexLabeled(
        cdsCondensedSequence::printSequence(cdsCondensedSequence::defaultConfigLabelled(), reordered));
}
