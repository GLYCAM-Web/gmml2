#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_PARSEDRESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_PARSEDRESIDUE_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/CodeUtils/casting.hpp"
#include <string>

namespace cdsCondensedSequence
{
    std::string getLink(cds::ResidueType type, const std::string& linkage);

    //	class ParsedResidue : public Abstract::absResidue , public glygraph::Node<ParsedResidue>
    class ParsedResidue : public cds::Residue
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ParsedResidue(const ParsedResidueComponents&);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::string GetName() const;
        std::string GetLinkageName() const;

        inline std::string GetInputString() const
        {
            return components.fullString;
        }

        inline std::string GetPreIsomerModifier() const
        {
            return components.preIsomerModifier;
        }

        inline std::string GetIsomer() const
        {
            return components.isomer;
        }

        inline std::string GetResidueName() const
        {
            return components.name;
        }

        inline std::string GetRingType() const
        {
            return components.ringType;
        }

        inline std::string GetRingShape() const
        {
            return components.ringShape;
        }

        inline std::string GetResidueModifier() const
        {
            return components.modifier;
        }

        inline std::string GetConfiguration() const
        {
            return components.configuration;
        }

        inline std::string GetLinkage() const
        {
            return components.linkage;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        inline std::string GetLink() const
        {
            return getLink(components.type, components.linkage);
        }

        std::vector<ParsedResidue*> GetChildren() const;
        std::vector<ParsedResidue*> GetParents() const;
        ParsedResidue* GetParent() const;
        std::string GetChildLinkagesForGlycamResidueNaming() const;

        //////////////////////////////////////////////////////////
        //                  OPERATOR OVERLOADING                //
        //////////////////////////////////////////////////////////
        bool operator==(const cds::Residue& rhs) const
        {
            return (this->GetLink() == codeUtils::erratic_cast<const ParsedResidue*>(&rhs)->GetLink());
        }

        bool operator!=(const cds::Residue& rhs) const
        {
            return (this->GetLink() != codeUtils::erratic_cast<const ParsedResidue*>(&rhs)->GetLink());
        }

        bool operator>(const cds::Residue& rhs) const
        {
            return (this->GetLink() > codeUtils::erratic_cast<const ParsedResidue*>(&rhs)->GetLink());
        }

        bool operator<(const cds::Residue& rhs) const
        {
            return (this->GetLink() < codeUtils::erratic_cast<const ParsedResidue*>(&rhs)->GetLink());
        }

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
        ParsedResidueComponents components;
    };
} // namespace cdsCondensedSequence
#endif
