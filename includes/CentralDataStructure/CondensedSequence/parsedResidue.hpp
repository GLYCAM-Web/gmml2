#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_PARSEDRESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_PARSEDRESIDUE_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/casting.hpp"
#include <string>

namespace cdsCondensedSequence
{
    using cds::ResidueType;

    //	class ParsedResidue : public Abstract::absResidue , public glygraph::Node<ParsedResidue>
    class ParsedResidue : public cds::Residue
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ParsedResidue(std::string residueString, ResidueType specifiedType = ResidueType::Undefined);
        ParsedResidue(std::string residueString, ParsedResidue* neighbor,
                      ResidueType specifiedType = ResidueType::Undefined);

        ~ParsedResidue()
        {} // std::cout << "ParsedResidue dtor for " << this->getName() << ", ";}

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::string GetName(const bool withLabels = false, const bool iupacConsensed = false) const;
        std::string GetIupacName() const;
        std::string GetLinkageName(const bool withLabels = false) const;

        inline std::string GetInputString() const
        {
            return fullResidueString_;
        }

        inline std::string GetIsomer() const
        {
            return isomer_;
        }

        inline std::string GetResidueName() const
        {
            return residueName_;
        }

        inline std::string GetRingType() const
        {
            return ringType_;
        }

        inline std::string GetRingShape() const
        {
            return ringShape_;
        }

        inline std::string GetResidueModifier() const
        {
            return residueModifier_;
        }

        inline std::string GetConfiguration() const
        {
            return configuration_;
        }

        inline std::string GetLinkage() const
        {
            return linkage_;
        }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void AddLinkage(ParsedResidue* otherRes);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string GetLink() const;
        std::vector<ParsedResidue*> GetChildren() const;
        std::vector<ParsedResidue*> GetParents() const;
        ParsedResidue* GetParent() const;
        std::string GetChildLinkagesForGlycamResidueNaming() const;
        std::string PrintToString() const;
        std::string GetMonosaccharideName() const;

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
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void ParseResidueStringIntoComponents(const std::string& residueString,
                                              ResidueType specifiedType = ResidueType::Undefined);
        void ExciseRingShapeFromModifier();

        //////////////////////////////////////////////////////////
        //                       MUTATORS                       //
        //////////////////////////////////////////////////////////
        inline void SetIsomer(std::string isomer)
        {
            isomer_ = isomer;
        }

        inline void SetResidueName(std::string name)
        {
            residueName_ = name;
        }

        inline void SetRingType(std::string type)
        {
            ringType_ = type;
        }

        inline void SetRingShape(std::string shape)
        {
            ringShape_ = shape;
        }

        inline void SetResidueModifier(std::string modifier)
        {
            residueModifier_ = modifier;
        }

        inline void SetConfiguration(std::string config)
        {
            configuration_ = config;
        }

        inline void SetLinkage(std::string label)
        {
            linkage_ = label;
        }

        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
        std::string fullResidueString_ = ""; // DManpNAca1-4, etc
        std::string isomer_            = ""; // D or L
        std::string residueName_       = ""; // Man, Neu, Ido etc
        std::string ringType_          = ""; // f or p
        std::string ringShape_         = ""; // 2SO, 4C1, 1C4 etc
        std::string residueModifier_   = ""; // NAc, Gc, A (IdoA) etc
        std::string configuration_     = ""; // a or b
        std::string linkage_           = ""; // 1-4, 2-6, 1- (when connected to OH) etc
    };
} // namespace cdsCondensedSequence
#endif
