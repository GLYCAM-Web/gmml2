#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_PARSEDRESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_PARSEDRESIDUE_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/casting.hpp"
#include <string>

namespace cdsCondensedSequence
{
    using cds::ResidueType;

    struct RingShapeAndModifier
    {
        std::string shape;
        std::string modifier;
    };

    struct ParsedResidueComponents
    {
        std::string fullString;
        cds::ResidueType type;
        std::string name;
        std::string linkage;
        std::string ringType;
        std::string configuration;
        std::string isomer;
        RingShapeAndModifier ringShapeAndModifier;
    };

    std::string getLink(const ParsedResidueComponents& components);

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
        std::string GetName(const bool withLabels = false, const bool iupacConsensed = false) const;
        std::string GetLinkageName(const bool withLabels = false) const;

        inline std::string GetInputString() const
        {
            return components.fullString;
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
            return components.ringShapeAndModifier.shape;
        }

        inline std::string GetResidueModifier() const
        {
            return components.ringShapeAndModifier.modifier;
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
            return getLink(components);
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
