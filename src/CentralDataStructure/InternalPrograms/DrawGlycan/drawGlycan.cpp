#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequencePrinter.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <string>

void CondensedSequence::drawGlycan(cdsCondensedSequence::GraphVizDotConfig configs, std::string condensedSequence)
{
    cdsCondensedSequence::printGraphViz(configs, cdsCondensedSequence::parseSequence(condensedSequence));
}
