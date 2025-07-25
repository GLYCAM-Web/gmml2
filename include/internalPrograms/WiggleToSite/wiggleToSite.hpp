#ifndef INCLUDE_INTERNALPROGRAMS_WIGGLETOSITE_WIGGLETOSITE_HPP
#define INCLUDE_INTERNALPROGRAMS_WIGGLETOSITE_WIGGLETOSITE_HPP

#include "include/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/internalPrograms/WiggleToSite/inputs.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/readers/Pdb/pdbFile.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/sequence/carbohydrate.hpp"
#include "include/util/containerTypes.hpp"

namespace gmml
{
    class WiggleToSite
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        WiggleToSite(const ParameterManager& parameterManager, WiggleToSiteInputs inputStruct);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        int minimizeDistance(
            const util::SparseVector<double>& elementRadii,
            const DihedralAngleDataTable& metadataTable,
            int persistCycles = 25,
            bool useMonteCarlo = true,
            int structureCount = 0);

      private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        Carbohydrate& getCarbohydrate() { return carbohydrate_; }

        std::vector<ResidueLinkage>& getWiggleLinkages() { return wiggleLinkages_; }

        pdb::PdbFile& getSubstrate() { return substrate_; }

        std::vector<Atom*> getAtomsToAvoid() { return atomsToAvoid_; }

        double getCurrentOverlap() { return currentOverlap_; }

        double getCurrentDistance() { return currentDistance_; }

        //////////////////////////////////////////////////////////
        //                  PRIVATE FUNCTIONS                   //
        //////////////////////////////////////////////////////////
        void setCurrentOverlap(double a) { currentOverlap_ = a; }

        void setCurrentDistance(double d) { currentDistance_ = d; }

        void superimpose(
            std::vector<Coordinate>& carbohydrateCoordinates,
            const Residue* superimpositionTarget,
            Residue* superimposeMe);
        std::vector<ResidueLinkage>& determineWiggleLinkages(
            const DihedralAngleDataTable& metadataTable, Residue* startResidue, Residue* endResidue);

        double calculateDistance();
        bool acceptOverlaps(const util::SparseVector<double>& elementRadii);
        bool acceptDistance(bool useMonteCarlo, double acceptance);
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        pdb::PdbFile substrate_;
        Carbohydrate carbohydrate_;
        std::vector<ResidueLinkage> wiggleLinkages_;
        std::vector<Atom*> atomsToAvoid_;
        std::vector<Atom*> wiggleMeAtoms_;
        std::vector<Atom*> wiggleTargetAtoms_;
        double currentOverlap_ = 0.0;
        double currentDistance_;
    };
} // namespace gmml

#endif
