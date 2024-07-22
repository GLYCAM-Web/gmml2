#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <array>
#include <vector>

namespace cds
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

    struct AngleOverlap
    {
        cds::Overlap overlaps;
        AngleWithMetadata angle;
    };

    struct DihedralRotationData
    {
        std::vector<Sphere> coordinates;
        std::vector<Sphere> boundingSpheres;
        std::vector<std::pair<size_t, size_t>> residueAtoms;
        std::vector<bool> firstResidueCoordinateMoving;
    };

    std::array<DihedralRotationData, 2> dihedralRotationInputData(RotatableDihedral& dihedral,
                                                                  const std::array<std::vector<Residue*>, 2>& residues);
    AngleOverlap wiggleWithinRangesDistanceCheck(RotatableDihedral& dihedral, std::vector<cds::Atom*>& overlapAtomSet1,
                                                 std::vector<cds::Atom*>& overlapAtomSet2,
                                                 const DihedralAngleData& metadata, std::vector<double> angles);
    AngleOverlap wiggleWithinCurrentRotamer(RotatableDihedral& dihedral, std::vector<cds::Atom*>& overlapAtomSet1,
                                            std::vector<cds::Atom*>& overlapAtomSet2, int angleIncrement);
    AngleOverlap wiggleUsingRotamers(const DihedralCoordinates coordinates, const DihedralAngleDataVector& rotamers,
                                     int angleIncrement, const std::array<DihedralRotationData, 2>& input);
    void simpleWiggleCurrentRotamers(std::vector<RotatableDihedral>& dihedrals,
                                     std::vector<cds::Atom*>& overlapAtomSet1, std::vector<cds::Atom*>& overlapAtomSet2,
                                     int angleIncrement);
    void simpleWiggleCurrentRotamers(std::vector<RotatableDihedral>& dihedrals,
                                     const std::array<std::vector<cds::Residue*>, 2>& residues, int angleIncrement);
} // namespace cds
#endif
