#ifndef INCLUDE_SEQUENCE_PARSEDRESIDUE_HPP
#define INCLUDE_SEQUENCE_PARSEDRESIDUE_HPP

#include "include/CentralDataStructure/residue.hpp"
#include "include/sequence/sequenceTypes.hpp"
#include "include/sequence/sequenceUtil.hpp"

#include <string>

namespace gmml
{
    //	class ParsedResidue : public Abstract::absResidue , public glygraph::Node<ParsedResidue>
    class ParsedResidue : public Residue
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ParsedResidue(const sequence::ParsedResidueComponents&);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::string GetName() const;
        std::string GetLinkageName() const;

        inline std::string GetInputString() const { return components.fullString; }

        inline std::string GetPreIsomerModifier() const { return components.preIsomerModifier; }

        inline std::string GetIsomer() const { return components.isomer; }

        inline std::string GetResidueName() const { return components.name; }

        inline std::string GetRingType() const { return components.ringType; }

        inline std::string GetRingShape() const { return components.ringShape; }

        inline std::string GetResidueModifier() const { return components.modifier; }

        inline std::string GetConfiguration() const { return components.configuration; }

        inline std::string GetLinkage() const { return components.linkage; }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        inline std::string GetLink() const { return sequence::mainLinkage(components.type, components.linkage); }

        std::vector<ParsedResidue*> GetChildren() const;
        std::vector<ParsedResidue*> GetParents() const;
        ParsedResidue* GetParent() const;

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
        sequence::ParsedResidueComponents components;
    };
} // namespace gmml

#endif
