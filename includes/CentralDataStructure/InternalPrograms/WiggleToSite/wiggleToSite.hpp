#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_WIGGLETOSITE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_WIGGLETOSITE_HPP

#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/inputs.hpp"
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/containerTypes.hpp"

namespace gmmlPrograms
{
    using cds::Residue;

    class WiggleToSite
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        WiggleToSite(const cdsParameters::ParameterManager& parameterManager, WiggleToSiteInputs inputStruct);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        int minimizeDistance(const codeUtils::SparseVector<double>& elementRadii,
                             const GlycamMetadata::DihedralAngleDataTable& metadataTable, int persistCycles = 25,
                             bool useMonteCarlo = true, int structureCount = 0);

      private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        cdsCondensedSequence::Carbohydrate& getCarbohydrate()
        {
            return carbohydrate_;
        }

        std::vector<cds::ResidueLinkage>& getWiggleLinkages()
        {
            return wiggleLinkages_;
        }

        pdb::PdbFile& getSubstrate()
        {
            return substrate_;
        }

        std::vector<cds::Atom*> getAtomsToAvoid()
        {
            return atomsToAvoid_;
        }

        cds::Overlap getCurrentOverlapCount()
        {
            return currentOverlapCount_;
        }

        double getCurrentDistance()
        {
            return currentDistance_;
        }

        //////////////////////////////////////////////////////////
        //                  PRIVATE FUNCTIONS                   //
        //////////////////////////////////////////////////////////
        void setCurrentOverlapCount(cds::Overlap i)
        {
            currentOverlapCount_ = i;
        }

        void setCurrentDistance(double d)
        {
            currentDistance_ = d;
        }

        void superimpose(std::vector<cds::CoordinateReference>& carbohydrateCoordinates,
                         const Residue* superimpositionTarget, Residue* superimposeMe);
        std::vector<cds::ResidueLinkage>&
        determineWiggleLinkages(const GlycamMetadata::DihedralAngleDataTable& metadataTable, Residue* startResidue,
                                Residue* endResidue);

        double calculateDistance();
        bool acceptOverlaps(const codeUtils::SparseVector<double>& elementRadii);
        bool acceptDistance(bool useMonteCarlo, double acceptance);
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        pdb::PdbFile substrate_;
        cdsCondensedSequence::Carbohydrate carbohydrate_;
        std::vector<cds::ResidueLinkage> wiggleLinkages_;
        std::vector<cds::Atom*> atomsToAvoid_;
        std::vector<cds::CoordinateReference> wiggleMeCoordinates_;
        std::vector<cds::CoordinateReference> wiggleTargetCoordinates_;
        cds::Overlap currentOverlapCount_ = cds::Overlap {0, 0.0};
        double currentDistance_;
    };
} // namespace gmmlPrograms
#endif
