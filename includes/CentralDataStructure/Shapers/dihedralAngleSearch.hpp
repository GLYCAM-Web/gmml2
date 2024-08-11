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

    struct AngleSearchPreference
    {
        std::vector<double> angles;
        std::vector<size_t> metadataOrder;
    };

    size_t bestOverlapResultIndex(const std::vector<AngleOverlap>& results);
    AngleOverlap bestOverlapResult(const std::vector<AngleOverlap>& results);
    std::array<DihedralRotationData, 2> dihedralRotationInputData(RotatableDihedral& dihedral,
                                                                  const std::array<std::vector<Residue*>, 2>& residues);
    AngleOverlap wiggleWithinRangesDistanceCheck(RotatableDihedral& dihedral, const AngleSearchPreference& preference,
                                                 std::vector<cds::Atom*>& overlapAtomSet1,
                                                 std::vector<cds::Atom*>& overlapAtomSet2, size_t metadataIndex,
                                                 std::vector<double> angles);
    AngleOverlap wiggleWithinCurrentRotamer(RotatableDihedral& dihedral, const DihedralAngleDataVector& metadataVector,
                                            const AngleSearchPreference& preference,
                                            std::vector<cds::Atom*>& overlapAtomSet1,
                                            std::vector<cds::Atom*>& overlapAtomSet2, double angleIncrement);
    AngleOverlap wiggleUsingRotamers(const DihedralCoordinates coordinates, const std::vector<size_t>& indices,
                                     const DihedralAngleDataVector& rotamers, const AngleSearchPreference& preference,
                                     double angleIncrement, const std::array<DihedralRotationData, 2>& input);
    void simpleWiggleCurrentRotamers(std::vector<RotatableDihedral>& dihedrals,
                                     const std::vector<DihedralAngleDataVector>& metadata,
                                     const std::vector<AngleSearchPreference>& preferences,
                                     std::vector<cds::Atom*>& overlapAtomSet1, std::vector<cds::Atom*>& overlapAtomSet2,
                                     double angleIncrement);
    void simpleWiggleCurrentRotamers(std::vector<RotatableDihedral>& dihedrals,
                                     const std::vector<DihedralAngleDataVector>& metadata,
                                     const std::vector<AngleSearchPreference>& preference,
                                     const std::array<std::vector<cds::Residue*>, 2>& residues, double angleIncrement);
    std::vector<AngleSearchPreference> angleSearchPreference(const ResidueLinkageShapePreference& preference);
    std::vector<std::vector<AngleSearchPreference>>
    angleSearchPreference(const std::vector<ResidueLinkageShapePreference>& preferences);

    // wiggle with:
    // shape preference angles
    // allowed min/max angles
    // metadata order (all for permutation, one for conformer)
    // within current rotamer overrides shape preference metadata order

} // namespace cds
#endif
