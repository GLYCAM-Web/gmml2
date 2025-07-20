#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearchTypes.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <string>
#include <variant>
#include <vector>

namespace cdsCondensedSequence
{
    class Carbohydrate : public cds::Molecule
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        Carbohydrate(
            const cdsParameters::ParameterManager& parameterManager,
            const codeUtils::SparseVector<double>& elementRadii,
            const SequenceData& sequence);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline std::string GetInputSequenceString() const { return inputSequenceString_; }

        inline std::vector<cds::ResidueLinkage>& GetGlycosidicLinkages() { return glycosidicLinkages_; }

        inline unsigned long int GetResidueCount() const { return this->getResidues().size(); }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void deleteResidue(cds::Residue* byeBye);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void Generate3DStructureFiles(
            const std::string& fileOutputDirectory,
            const std::string& outputFileNaming,
            const std::vector<std::string>& headerLines);
        void ResolveOverlaps(
            const codeUtils::SparseVector<double>& elementRadii,
            const GlycamMetadata::DihedralAngleDataTable& metadataTable,
            const cds::AngleSearchSettings& searchSettings);
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
