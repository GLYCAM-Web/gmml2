#include "includes/CentralDataStructure/CondensedSequence/sequenceParser2.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"

using cdsCondensedSequence::Sequence;

Sequence::Sequence(std::string glycamCondensedString) : cds::Molecule()
{
    parseSequence(this, glycamCondensedString);
    this->setName("CONDENSEDSEQUENCE");
    reorderSequence(this);
}

void Sequence::LabelSequence()
{
    labelSequence(this->getResidues());
}

std::string Sequence::Print(const bool withLabels, const bool iupac) const
{
    return printSequence(this->getResidues(), withLabels, iupac);
}

std::string Sequence::PrintGraphViz(GraphVizDotConfig& configs)
{
    return printGraphViz(configs, this->getResidues());
}
