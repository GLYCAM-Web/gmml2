#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include <variant>
#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    constexpr auto defaultSearchAngles =
        [](const GlycamMetadata::DihedralAngleData& metadata, double preference, double deviation)
    {
        double increment = 1.0;
        auto angle       = metadata.angle_deviation;
        if (std::holds_alternative<GlycamMetadata::AngleLimit>(angle))
        {
            auto dev = std::get<GlycamMetadata::AngleLimit>(angle);
            return cds::evenlySpacedAngles(preference, dev.lowerDeviationLimit, dev.upperDeviationLimit, increment);
        }
        else if (std::holds_alternative<GlycamMetadata::AngleStd>(angle))
        {
            auto dev = std::get<GlycamMetadata::AngleStd>(angle);
            return cds::evenlySpacedAngles(preference, deviation * dev.lowerDeviationStd,
                                           deviation * dev.upperDeviationStd, increment);
        }
        else
        {
            throw std::runtime_error("unknown dihedral angle type");
        }
    };
    const cds::AngleSearchSettings defaultSearchSettings = {1.0, defaultSearchAngles};

    class Carbohydrate : public cds::Molecule
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        Carbohydrate(std::string inputSequence);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline std::string GetInputSequenceString() const
        {
            return inputSequenceString_;
        }

        inline std::vector<cds::ResidueLinkage>& GetGlycosidicLinkages()
        {
            return glycosidicLinkages_;
        }

        inline unsigned long int GetResidueCount() const
        {
            return this->getResidues().size();
        }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void deleteResidue(cds::Residue* byeBye);
        void replaceAglycone(cds::Residue* newAglycone);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void Generate3DStructureFiles(std::string fileOutputDirectory = "unspecified",
                                      std::string outputFileNaming    = "structure");
        void ResolveOverlaps(const cds::AngleSearchSettings& searchSettings);
        unsigned long int CountShapes(bool likelyShapesOnly = false) const;
        std::string GetNumberOfShapes(
            bool likelyShapesOnly = false) const; // This one is for gems. ToDo try to deprecate and use CountShapes.
        cds::Residue* GetReducingResidue();
        cds::Residue* GetAglycone();
        cds::Atom* GetAnomericAtom();

      private:
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void ApplyDeoxy(ParsedResidue* deoxyResidue);
        void DerivativeChargeAdjustment(ParsedResidue* parsedResidue);
        void ConnectAndSetGeometry(cds::Residue* parentResidue, cds::Residue* childResidue,
                                   const cds::AngleSearchSettings& searchSettings);
        std::vector<std::string> GetGlycamNamesOfResidues() const;
        std::string GetGlycamResidueName(ParsedResidue* residue) const;
        void DepthFirstSetConnectivityAndGeometry(cds::Residue* currentParent,
                                                  const cds::AngleSearchSettings& searchSettings);
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        std::string inputSequenceString_;
        std::vector<cds::ResidueLinkage> glycosidicLinkages_;
    };
} // namespace cdsCondensedSequence
#endif
