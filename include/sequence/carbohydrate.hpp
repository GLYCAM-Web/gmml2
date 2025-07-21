#ifndef INCLUDES_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP
#define INCLUDES_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP

#include "include/CentralDataStructure/Shapers/dihedralAngleSearchTypes.hpp"
#include "include/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/sequence/sequenceTypes.hpp"
#include "include/util/containerTypes.hpp"

#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    class Carbohydrate : public Molecule
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        Carbohydrate(
            const ParameterManager& parameterManager,
            const util::SparseVector<double>& elementRadii,
            const sequence::SequenceData& sequence);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline std::string GetInputSequenceString() const { return inputSequenceString_; }

        inline std::vector<ResidueLinkage>& GetGlycosidicLinkages() { return glycosidicLinkages_; }

        inline unsigned long int GetResidueCount() const { return this->getResidues().size(); }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void deleteResidue(Residue* byeBye);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void Generate3DStructureFiles(
            const std::string& fileOutputDirectory,
            const std::string& outputFileNaming,
            const std::vector<std::string>& headerLines);
        void ResolveOverlaps(
            const util::SparseVector<double>& elementRadii,
            const DihedralAngleDataTable& metadataTable,
            const AngleSearchSettings& searchSettings);
        unsigned long int CountShapes(bool likelyShapesOnly = false) const;
        std::string GetNumberOfShapes(
            bool likelyShapesOnly = false) const; // This one is for gems. ToDo try to deprecate and use CountShapes.

      private:
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        std::string inputSequenceString_;
        std::vector<ResidueLinkage> glycosidicLinkages_;
    };
} // namespace gmml

#endif
