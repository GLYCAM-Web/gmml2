#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <memory>
#include <string>
#include <vector>

using CondensedSequence::DrawGlycan;

DrawGlycan::DrawGlycan(cdsCondensedSequence::GraphVizDotConfig configs, std::string condensedSequence)
{
    std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>> residues;
    cdsCondensedSequence::parseSequence(residues, condensedSequence);
    cdsCondensedSequence::printGraphViz(configs, codeUtils::pointerToUniqueVector(residues));
}
