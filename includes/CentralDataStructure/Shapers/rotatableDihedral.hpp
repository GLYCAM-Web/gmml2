#ifndef GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_ROTATABLE_DIHEDRAL_HPP
#define GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_ROTATABLE_DIHEDRAL_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/CentralDataStructure/Overlaps/cdsOverlaps.hpp"

#include <random>

using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;
// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32
    rng(seed_source); // Eclipse complains about ambiguity, and yet it compiles... // @suppress("Ambiguous problem")

// This class stores the four atoms that define a dihedral angle, the atoms that move when it is rotated
// and, if moved, the previous dihedral angle, which allows me to reset easily.
namespace cds
{
    struct AngleOverlap
    {
        double angle;
        unsigned int overlaps;
    };

    struct AngleWithMetadata
    {
        double angle;
        const DihedralAngleData* metadata;
    };

    class RotatableDihedral
    {
      public:
        RotatableDihedral(const std::array<Atom*, 4>& atoms, const std::vector<DihedralAngleData>& metadata)
            : atoms_(atoms), assigned_metadata_(metadata) {};

        ////////////////////////////////////////////
        //////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline const DihedralAngleDataVector& GetMetadata() const
        {
            return assigned_metadata_;
        }

        DihedralAngleDataVector GetLikelyMetadata() const;
        int GetNumberOfRotamers(bool likelyShapesOnly = false) const;
        std::string GetName() const;
        double CalculateDihedralAngle() const;
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void DetermineAtomsThatMove();     // Based on connectivities, this figures out which atoms will move when the
                                           // dihedral is rotated.
        void SetDihedralAngleToPrevious(); // Sets the dihedral to previous dihedral angle

        void SetRandomAngleEntryUsingMetadata();
        void SetSpecificAngleEntryUsingMetadata(bool useRanges, long unsigned int angleEntryNumber);
        bool SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
        void WiggleWithinCurrentRotamer(std::vector<cds::Atom*>& overlapAtomSet1,
                                        std::vector<cds::Atom*>& overlapAtomSet2, int angleIncrement);
        void WiggleWithinCurrentRotamer(std::vector<cds::Residue*>& overlapResidueSet1,
                                        std::vector<cds::Residue*>& overlapResidueSet2, int angleIncrement);
        void WiggleUsingAllRotamers(std::vector<cds::Residue*>& overlapAtomSet1,
                                    std::vector<cds::Residue*>& overlapAtomSet2, int angleIncrement);
        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        std::string Print() const;

      private:
        //////////////////////////////////////////////////////////
        //                  PRIVATE ACCESSORS                   //
        //////////////////////////////////////////////////////////
        inline AngleWithMetadata GetPreviousState() const
        {
            return previousState_;
        }

        //////////////////////////////////////////////////////////
        //                  PRIVATE MUTATORS                    //
        //////////////////////////////////////////////////////////
        void SetAtomsThatMove(std::vector<cds::Atom*>& atoms);
        void SetDihedralAngle(
            double dihedral_angle,
            const DihedralAngleData* metadata); // Sets the dihedral angle by rotating the bond between atom2 and
        double RandomizeDihedralAngleWithinMetadataRange(
            const DihedralAngleData* entry); // Randomly sets dihedral angle to a value allowed by the metadata

        //////////////////////////////////////////////////////////
        //                  PRIVATE FUNCTIONS                   //
        //////////////////////////////////////////////////////////

        std::array<Coordinate*, 4> dihedralCoordinates() const;

        inline void RecordPreviousState(double d, const DihedralAngleData* metadata)
        {
            previousState_ = {d, metadata};
        }

        inline void SetWasEverRotated(bool b)
        {
            wasEverRotated_ = b;
        }

        inline void SetCurrentMetaData(const DihedralAngleData* d)
        {
            currentMetadata_ = d;
        }

        inline const DihedralAngleData* GetCurrentMetaData()
        {
            return currentMetadata_;
        }

        std::vector<cds::AngleOverlap> WiggleWithinRangesDistanceCheck(std::vector<cds::Atom*>& overlapAtomSet1,
                                                                       std::vector<cds::Atom*>& overlapAtomSet2,
                                                                       const DihedralAngleData* metadata,
                                                                       std::vector<double> angles);
        std::vector<cds::AngleOverlap> WiggleWithinRangesDistanceCheck(cds::ResidueAtomOverlapInputPair& overlapInput,
                                                                       const DihedralAngleData* metadata,
                                                                       std::vector<double> angles);

        inline std::vector<cds::Coordinate*>& GetCoordinatesThatMove()
        {
            return coordinatesThatMove_;
        }

        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
        std::array<cds::Atom*, 4> atoms_;
        // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond
        // is rotated.
        //    std::vector<cds::Atom*> atoms_that_move_;
        //    std::vector<cds::Atom*> extra_atoms_that_move_;
        std::vector<cds::Coordinate*> coordinatesThatMove_;
        AngleWithMetadata previousState_ = {constants::dNotSet,
                                            nullptr}; // I often want to reset a dihedral angle after rotating
                                                      // it, so recording the previous angle makes this easy.
        DihedralAngleDataVector assigned_metadata_;
        const DihedralAngleData* currentMetadata_ = nullptr;
        bool wasEverRotated_                      = false; // Need this, as it might add a H atom for psi
    };

} // namespace cds
#endif
