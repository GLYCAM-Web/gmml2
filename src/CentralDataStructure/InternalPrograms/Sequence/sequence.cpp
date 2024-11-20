#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <memory>
#include <vector>

using CondensedSequence::Sequence;

Sequence::Sequence(std::string condensedSequence)
{
    std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>> residuePtrs;
    cdsCondensedSequence::parseSequence(residuePtrs, condensedSequence);
    this->setInputSequence(condensedSequence);
    this->setInterpretedSequence(
        cdsCondensedSequence::printSequence(codeUtils::pointerToUniqueVector(residuePtrs), false));
    this->setIndexOrdered(cdsCondensedSequence::reorderSequence(residuePtrs));
    std::vector<cdsCondensedSequence::ParsedResidue*> residues = codeUtils::pointerToUniqueVector(residuePtrs);
    cdsCondensedSequence::labelSequence(residues);
    this->setIndexLabeled(cdsCondensedSequence::printSequence(residues, true));
}
