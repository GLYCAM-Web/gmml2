#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"

#include <string>

using CondensedSequence::DrawGlycan;

DrawGlycan::DrawGlycan(cdsCondensedSequence::GraphVizDotConfig configs, std::string condensedSequence)
{
    cds::Molecule molecule;
    cdsCondensedSequence::parseSequence(&molecule, condensedSequence);
    cdsCondensedSequence::printGraphViz(configs, molecule.getResidues());
}
