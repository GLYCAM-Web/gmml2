#include "include/sequence/sequence.hpp"

#include "include/sequence/sequenceManipulation.hpp"
#include "include/sequence/sequenceParser.hpp"
#include "include/sequence/sequencePrinter.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <vector>

namespace gmml
{
    namespace sequence
    {
        Sequence::Sequence(std::string condensedSequence)
        {
            AbstractSequence abstract = parseSequence(condensedSequence);
            SequenceData sequence = instantiate(abstract);
            this->setInterpretedSequence(printSequence(defaultConfigUnsorted(), sequence));
            SequenceData reordered = sequence::reordered(sequence);
            this->setInputSequence(condensedSequence);
            this->setIndexOrdered(printSequence(defaultConfig(), reordered));
            this->setIndexLabeled(printSequence(defaultConfigLabelled(), reordered));
        }
    } // namespace sequence
} // namespace gmml
