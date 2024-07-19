#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATABLEDIHEDRAL_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATABLEDIHEDRAL_HPP

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

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
    struct AngleWithMetadata
    {
        double value;
        const DihedralAngleData* metadata;
    };

    struct AngleOverlap
    {
        cds::Overlap overlaps;
        AngleWithMetadata angle;
    };

    struct dihedralRotationData
    {
        std::vector<Sphere> coordinates;
        std::vector<Sphere> boundingSpheres;
        std::vector<std::pair<size_t, size_t>> residueAtoms;
        std::vector<bool> firstResidueCoordinateMoving;
    };

    std::array<dihedralRotationData, 2> dihedralRotationInputData(bool branched, const std::array<Atom*, 4>& dihedral,
                                                                  const std::vector<Coordinate*>& movingCoordinates,
                                                                  const std::array<std::vector<Residue*>, 2>& residues);

    class RotatableDihedral
    {
      public:
        RotatableDihedral(bool isBranchingLinkage, const std::array<Atom*, 4>& atoms,
                          const std::vector<DihedralAngleData>& metadata)
            : isBranchingLinkage_(isBranchingLinkage), atoms_(atoms), assigned_metadata_(metadata) {};

        ////////////////////////////////////////////
        //////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline const DihedralAngleDataVector& GetMetadata() const
        {
            return assigned_metadata_;
        }

        inline const DihedralAngleData* GetCurrentMetaData()
        {
            return currentMetadata_;
        }

        inline AngleWithMetadata GetPreviousState() const
        {
            return previousState_;
        }

        DihedralAngleDataVector GetLikelyMetadata() const;
        int GetNumberOfRotamers(bool likelyShapesOnly = false) const;
        std::string GetName() const;
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void DetermineAtomsThatMove(); // Based on connectivities, this figures out which atoms will move when the
                                       // dihedral is rotated.
        void
        SetDihedralAngle(AngleWithMetadata target); // Sets the dihedral angle by rotating the bond between atom2 and

        AngleWithMetadata RandomAngleEntryUsingMetadata();
        AngleWithMetadata SpecificAngleEntryUsingMetadata(bool useRanges, long unsigned int angleEntryNumber);
        bool SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
        AngleOverlap WiggleWithinCurrentRotamer(std::vector<cds::Atom*>& overlapAtomSet1,
                                                std::vector<cds::Atom*>& overlapAtomSet2, int angleIncrement);
        AngleOverlap WiggleUsingRotamers(const DihedralAngleDataVector& rotamers, int angleIncrement,
                                         const std::array<std::vector<cds::Residue*>, 2>& residues);
        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        std::string Print() const;

      private:
        //////////////////////////////////////////////////////////
        //                  PRIVATE FUNCTIONS                   //
        //////////////////////////////////////////////////////////

        DihedralCoordinates dihedralCoordinates() const;

        inline void RecordPreviousState(AngleWithMetadata angle)
        {
            previousState_ = angle;
        }

        inline void SetCurrentMetaData(const DihedralAngleData* d)
        {
            currentMetadata_ = d;
        }

        cds::AngleOverlap WiggleWithinRangesDistanceCheck(std::vector<cds::Atom*>& overlapAtomSet1,
                                                          std::vector<cds::Atom*>& overlapAtomSet2,
                                                          const DihedralAngleData* metadata,
                                                          std::vector<double> angles);

        inline std::vector<cds::Coordinate*>& GetCoordinatesThatMove()
        {
            return coordinatesThatMove_;
        }

        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        bool isBranchingLinkage_;
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
    };

} // namespace cds
#endif
