#include "includes/CentralDataStructure/CondensedSequence/sequenceParser2.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CodeUtils/containers.hpp"

using cdsCondensedSequence::Sequence;

Sequence::Sequence(std::string glycamCondensedString)
{
    parseSequence(residues, glycamCondensedString);
    reorderSequence(residues);
}

void Sequence::LabelSequence()
{
    labelSequence(codeUtils::pointerToUniqueVector(residues));
}

std::string Sequence::Print(const bool withLabels)
{
    return printSequence(codeUtils::pointerToUniqueVector(residues), withLabels);
}

std::string Sequence::PrintIupac()
{
    return printSequence(codeUtils::pointerToUniqueVector(residues), false, true);
}

std::string Sequence::PrintGraphViz(GraphVizDotConfig& configs)
{
    return printGraphViz(configs, codeUtils::pointerToUniqueVector(residues));
}
