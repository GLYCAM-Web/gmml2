#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"

#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    class SequenceParser : public cds::Molecule
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        SequenceParser(std::string sequence);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::string Print() const;
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
      private:
        SequenceParser() {};

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline bool DerivativesExist()
        {
            return savedDerivatives_.size();
        }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        std::vector<std::string> ExtractDerivatives();
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        bool CheckSequenceSanity(std::string sequence);
        void ParseLabelledInput(std::string inString);
        std::string parseRepeatingUnits(const std::string inputSequence);
        unsigned int seekRepeatStart(const std::string& inputSequence, unsigned int startPosition);
        bool ParseCondensedSequence(std::string inString);
        void RecurveParseAlt(size_t& currentIndex, std::string sequence, ParsedResidue* parent);
        ParsedResidue* SaveResidue(const size_t windowStart, const size_t windowEnd, const std::string sequence,
                                   ParsedResidue* parent);

        inline void SaveDerivative(std::string derivative)
        {
            savedDerivatives_.push_back(derivative);
        }

        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        std::vector<std::string> savedDerivatives_;
    };
} // namespace cdsCondensedSequence
#endif
