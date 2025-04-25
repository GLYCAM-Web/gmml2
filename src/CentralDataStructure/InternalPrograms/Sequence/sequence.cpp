#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>

CondensedSequence::Sequence::Sequence(std::string condensedSequence)
{
    cdsCondensedSequence::SequenceData sequence = cdsCondensedSequence::parseSequence(condensedSequence);
    this->setInterpretedSequence(
        cdsCondensedSequence::printSequence(sequence, cdsCondensedSequence::defaultConfigUnsorted()));
    cdsCondensedSequence::SequenceData reordered =
        cdsCondensedSequence::reordered(cdsCondensedSequence::parseSequence(condensedSequence));
    this->setInputSequence(condensedSequence);
    this->setIndexOrdered(cdsCondensedSequence::printSequence(reordered, cdsCondensedSequence::defaultConfig()));
    this->setIndexLabeled(
        cdsCondensedSequence::printSequence(reordered, cdsCondensedSequence::defaultConfigLabelled()));
}
