#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    constexpr auto defaultSearchAngles =
        [](const GlycamMetadata::DihedralAngleData& metadata, double preference, double deviation)
    {
        double increment = 1.0;
        return cds::evenlySpacedAngles(preference, deviation, increment, metadata);
    };
    const cds::AngleSearchSettings defaultSearchSettings = {2.0, defaultSearchAngles};

    class Carbohydrate : public SequenceManipulator
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        Carbohydrate(std::string inputSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-"
                                                 "4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH");

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
        void SetDefaultShapeUsingMetadata();
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
