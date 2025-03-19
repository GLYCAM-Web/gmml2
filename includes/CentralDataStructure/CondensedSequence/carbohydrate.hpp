#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP

#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include <functional>
#include <variant>
#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    constexpr auto defaultSearchAngles =
        [](const GlycamMetadata::DihedralAngleData& metadata, double preference, double deviation)
    {
        double increment = 1.0;
        std::function<std::vector<double>(const GlycamMetadata::AngleLimit&)> onLimit =
            [&](const GlycamMetadata::AngleLimit& dev)
        {
            return cds::evenlySpacedAngles(preference, dev.lowerDeviationLimit, dev.upperDeviationLimit, increment);
        };
        std::function<std::vector<double>(const GlycamMetadata::AngleStd&)> onStd =
            [&](const GlycamMetadata::AngleStd& dev)
        {
            return cds::evenlySpacedAngles(preference, deviation * dev.lowerDeviationStd,
                                           deviation * dev.upperDeviationStd, increment);
        };
        return GlycamMetadata::onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
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
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void Generate3DStructureFiles(const std::string& fileOutputDirectory, const std::string& outputFileNaming,
                                      const std::vector<std::string>& headerLines);
        void ResolveOverlaps(const cds::AngleSearchSettings& searchSettings);
        unsigned long int CountShapes(bool likelyShapesOnly = false) const;
        std::string GetNumberOfShapes(
            bool likelyShapesOnly = false) const; // This one is for gems. ToDo try to deprecate and use CountShapes.

      private:
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        std::string inputSequenceString_;
        std::vector<cds::ResidueLinkage> glycosidicLinkages_;
    };
} // namespace cdsCondensedSequence
#endif
