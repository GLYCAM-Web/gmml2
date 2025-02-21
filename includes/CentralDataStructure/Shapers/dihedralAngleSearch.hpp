#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"

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

    struct DihedralRotationDataContainer
    {
        assembly::Graph graph;
        assembly::Bounds bounds;
        std::vector<bool> atomMoving;
        std::vector<bool> atomIncluded;
        std::vector<MolecularMetadata::Element> atomElements;
        std::vector<double> residueWeights;
        std::array<std::vector<size_t>, 2> residueIndices;
        std::vector<cds::BondedResidueOverlapInput> bonds;
    };

    struct DihedralRotationData
    {
        const assembly::Graph& graph;
        const assembly::Bounds& bounds;
        const std::vector<bool>& atomMoving;
        const std::vector<bool>& atomIncluded;
        const std::vector<MolecularMetadata::Element>& atomElements;
        const std::vector<double>& residueWeights;
        const std::array<std::vector<size_t>, 2>& residueIndices;
        const std::vector<cds::BondedResidueOverlapInput>& bonds;
    };

    struct AngleSearchPreference
    {
        double deviation;
        std::vector<double> angles;
        std::vector<size_t> metadataOrder;
    };

    typedef std::function<std::vector<double>(const DihedralAngleData&, double, double)> SearchAngles;

    struct AngleSearchSettings
    {
        double deviation;
        SearchAngles angles;
    };

    size_t bestOverlapResultIndex(const std::vector<AngleOverlap>& results);
    AngleOverlap bestOverlapResult(const std::vector<AngleOverlap>& results);
    DihedralRotationDataContainer dihedralRotationInputData(double overlapTolerance, RotatableDihedral& dihedral,
                                                            const GraphIndexData& indices,
                                                            const std::array<ResiduesWithOverlapWeight, 2>& residues);
    AngleOverlap wiggleUsingRotamers(const MolecularMetadata::PotentialTable& potential,
                                     cds::OverlapProperties overlapProperties, SearchAngles searchAngles,
                                     const DihedralCoordinates coordinates, const std::vector<size_t>& indices,
                                     const DihedralAngleDataVector& rotamers, const AngleSearchPreference& preference,
                                     const DihedralRotationData& input);
    void simpleWiggleCurrentRotamers(const MolecularMetadata::PotentialTable& potential,
                                     cds::OverlapProperties overlapProperties, SearchAngles searchAngles,
                                     std::vector<RotatableDihedral>& dihedrals,
                                     const std::vector<DihedralAngleDataVector>& metadata,
                                     const std::vector<AngleSearchPreference>& preference,
                                     const GraphIndexData& indices,
                                     const std::array<ResiduesWithOverlapWeight, 2>& residues);
    std::vector<double> evenlySpacedAngles(double preference, double lowerDeviation, double upperDeviation,
                                           double increment);
    std::vector<AngleSearchPreference> angleSearchPreference(double deviation,
                                                             const ResidueLinkageShapePreference& preference);
    std::vector<std::vector<AngleSearchPreference>>
    angleSearchPreference(double deviation, const std::vector<ResidueLinkageShapePreference>& preferences);
} // namespace cds
#endif
