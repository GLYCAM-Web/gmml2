#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"

#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulation.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequencePrinter.hpp"

#include <string>

void CondensedSequence::drawGlycan(cdsCondensedSequence::GraphVizDotConfig configs, std::string condensedSequence)
{
    cdsCondensedSequence::printGraphViz(
        configs, cdsCondensedSequence::instantiate(cdsCondensedSequence::parseSequence(condensedSequence)));
}
