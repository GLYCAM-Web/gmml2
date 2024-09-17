#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CodeUtils/logging.hpp"

using CondensedSequence::Sequence;

Sequence::Sequence(std::string condensedSequence)
{
    cds::Molecule molecule;
    cdsCondensedSequence::parseSequence(&molecule, condensedSequence);
    this->setInputSequence(condensedSequence);
    this->setInterpretedSequence(cdsCondensedSequence::printSequence(molecule.getResidues(), false));
    this->setIndexOrdered(cdsCondensedSequence::reorderSequence(&molecule));
    cdsCondensedSequence::labelSequence(molecule.getResidues());
    this->setIndexLabeled(cdsCondensedSequence::printSequence(molecule.getResidues(), true));
}
