#include "include/programs/swig/sequence.hpp"

#include "include/sequence/sequenceManipulation.hpp"
#include "include/sequence/sequenceParser.hpp"
#include "include/sequence/sequencePrinter.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <vector>

namespace gmml
{
    using namespace sequence;

    Sequence::Sequence(std::string condensedSequence)
    {
        AbstractSequence abstract = parseSequence(condensedSequence);
        SequenceData sequence = instantiate(abstract);
        interpretedSequence = printSequence(defaultConfigUnsorted(), sequence);
        SequenceData reordered = sequence::reordered(sequence);
        indexOrdered = printSequence(defaultConfig(), reordered);
        indexLabeled = printSequence(defaultConfigLabelled(), reordered);
    }
} // namespace gmml
