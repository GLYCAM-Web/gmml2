#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <array>
#include <vector>
#include <functional>

namespace cds
{
    using GlycamMetadata::DihedralAngleData;
    using GlycamMetadata::DihedralAngleDataVector;

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
        std::vector<double> residueWeights;
        std::vector<bool> firstResidueCoordinateMoving;
        std::vector<bool> firstResidueBondedAtoms;
    };

    struct AngleSearchPreference
    {
        std::vector<double> angles;
        std::vector<size_t> metadataOrder;
    };

    typedef std::function<std::vector<double>(const DihedralAngleData&)> SearchAngles;

    size_t bestOverlapResultIndex(const std::vector<AngleOverlap>& results);
    AngleOverlap bestOverlapResult(const std::vector<AngleOverlap>& results);
    std::array<DihedralRotationData, 2>
    dihedralRotationInputData(RotatableDihedral& dihedral, const std::array<ResiduesWithOverlapWeight, 2>& residues);
    AngleOverlap wiggleUsingRotamers(SearchAngles searchAngles, const DihedralCoordinates coordinates,
                                     const std::vector<size_t>& indices, const DihedralAngleDataVector& rotamers,
                                     const AngleSearchPreference& preference,
                                     const std::array<DihedralRotationData, 2>& input);
    void simpleWiggleCurrentRotamers(SearchAngles searchAngles, std::vector<RotatableDihedral>& dihedrals,
                                     const std::vector<DihedralAngleDataVector>& metadata,
                                     const std::vector<AngleSearchPreference>& preference,
                                     const std::array<ResiduesWithOverlapWeight, 2>& residues);
    std::vector<double> evenlySpacedAngles(double deviation, double increment, const DihedralAngleData& metadata);
    std::vector<AngleSearchPreference> angleSearchPreference(const ResidueLinkageShapePreference& preference);
    std::vector<std::vector<AngleSearchPreference>>
    angleSearchPreference(const std::vector<ResidueLinkageShapePreference>& preferences);
} // namespace cds
#endif
