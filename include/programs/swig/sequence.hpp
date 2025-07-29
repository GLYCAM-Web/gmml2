#ifndef INCLUDE_SEQUENCE_SEQUENCE_HPP
#define INCLUDE_SEQUENCE_SEQUENCE_HPP

// OG Feb 2022
// This class only exists so I can wrap it into gems without swig having to know about SequenceManipulator and the
// template Node class I don't know how to wrap that class, so I'm trying to keep it "on the other side of the line"
// when wrapping.

#include <string>

namespace gmml
{
    class Sequence
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTORS                   //
        //////////////////////////////////////////////////////////
        Sequence(std::string condensedSequence);

        inline std::string getInterpretedSequence() { return interpretedSequence; }

        inline std::string getIndexOrdered() { return indexOrdered; }

        inline std::string getIndexLabeled() { return indexLabeled; }

      private:
        std::string interpretedSequence = "";
        std::string indexOrdered = "";
        std::string indexLabeled = "";
    };
} // namespace gmml

#endif
