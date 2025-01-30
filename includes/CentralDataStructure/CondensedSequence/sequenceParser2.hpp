#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER2_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER2_HPP_

#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"

#include <memory>
#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    class Sequence
    {
      public:
        Sequence(std::string glycamCondensedString);
        void LabelSequence();
        std::string Print(const bool withLabels = false);
        std::string PrintIupac();
        std::string PrintGraphViz(GraphVizDotConfig& configs);

      private:
        std::string inputSequenceString_;
        std::vector<std::unique_ptr<ParsedResidue>> residues;
    };

} // namespace cdsCondensedSequence

#endif /* INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER2_HPP_ */
