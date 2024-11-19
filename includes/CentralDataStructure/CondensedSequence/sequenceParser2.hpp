#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER2_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER2_HPP_

#include <string>
#include "includes/CentralDataStructure/CondensedSequence/graphVizDotConfig.hpp"
#include "includes/CentralDataStructure/molecule.hpp"

namespace cdsCondensedSequence
{
    class Sequence : public cds::Molecule
    {
      public:
        Sequence(std::string glycamCondensedString);
        void LabelSequence();
        std::string Print(const bool withLabels = false) const;
        std::string PrintIupac() const;
        std::string PrintGraphViz(GraphVizDotConfig& configs);

      private:
        std::string inputSequenceString_;
    };

} // namespace cdsCondensedSequence

#endif /* INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER2_HPP_ */
