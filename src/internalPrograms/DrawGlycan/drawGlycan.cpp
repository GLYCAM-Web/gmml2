#include "include/internalPrograms/DrawGlycan/drawGlycan.hpp"

#include "include/sequence/sequenceManipulation.hpp"
#include "include/sequence/sequenceParser.hpp"
#include "include/sequence/sequencePrinter.hpp"

#include <string>

namespace gmml
{
    void drawGlycan(sequence::GraphVizDotConfig configs, std::string condensedSequence)
    {
        printGraphViz(configs, sequence::instantiate(sequence::parseSequence(condensedSequence)));
    }
} // namespace gmml
