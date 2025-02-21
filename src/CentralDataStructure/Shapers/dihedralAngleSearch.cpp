#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"

#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <cmath>
#include <numeric>
#include <vector>

using GlycamMetadata::DihedralAngleData;

namespace
{
    std::vector<double> evenlySpaced(double lower, double upper, double approximateIncrement)
    {
        double range = upper - lower;
        int steps    = std::ceil(std::abs(range) / approximateIncrement);
        if (steps == 0)
        {
            // range == 0, hence lower == upper
            return {};
        }
        else
        {
            double increment = range / steps;
            std::vector<double> result;
            result.reserve(steps);
            for (int k = 0; k < steps; k++)
            {
                result.push_back(lower + k * increment);
            }
            return result;
        }
    }

    std::vector<size_t> toResidueIndices(const std::vector<cds::Residue*>& residues,
                                         const std::vector<cds::Residue*>& reindex)
    {
        std::vector<size_t> result;
        result.reserve(reindex.size());
        for (auto res : reindex)
        {
            result.push_back(codeUtils::indexOf(residues, res));
        }
        return result;
    }

    std::array<std::vector<size_t>, 2> branchedResidueSets(const std::vector<bool>& residueMoving,
                                                           const std::vector<size_t>& setA,
                                                           const std::vector<size_t>& setB)
    {
        std::vector<size_t> newSetA = setB;
        std::vector<size_t> newSetB = {setA[0]};
        for (size_t n = 1; n < setA.size(); n++)
        {
            size_t index = setA[n];
            if (residueMoving[index])
            {
                newSetB.push_back(index);
            }
            else
            {
                newSetA.push_back(index);
            }
        }
        return {newSetA, newSetB};
    }

    std::vector<bool> toAtomMoving(const std::vector<cds::Atom*>& atoms, const std::vector<cds::Atom*>& movingAtoms)
    {
        std::vector<bool> atomMoving;
        atomMoving.reserve(atoms.size());
        for (auto& atom : atoms)
        {
            atomMoving.push_back(codeUtils::contains(movingAtoms, atom));
        }
        return atomMoving;
    }

    void applyMatrix(const assembly::Graph& graph, const assembly::Bounds& initial,
                     const std::vector<size_t>& movingAtoms, const cds::RotationMatrix& matrix,
                     assembly::Bounds& bounds)
    {
        for (size_t n : movingAtoms)
        {
            bounds.atoms[n].center = matrix * initial.atoms[n].center;
        }
        assembly::updateBoundsContainingAtoms(graph, bounds, movingAtoms);
    };

    cds::AngleOverlap WiggleAnglesOverlaps(const MolecularMetadata::PotentialTable& potential,
                                           cds::OverlapProperties overlapProperties,
                                           const cds::DihedralCoordinates dihedral, size_t metadataIndex,
                                           double anglePreference, std::vector<double> angles,
                                           const cds::DihedralRotationData& input)
    {
        const std::vector<size_t>& fixedResidueIndices  = input.residueIndices[0];
        const std::vector<size_t>& movingResidueIndices = input.residueIndices[1];
        std::vector<size_t> movingAtoms                 = codeUtils::boolsToIndices(input.atomMoving);
        assembly::Bounds bounds                         = input.bounds;
        std::vector<cds::AngleOverlap> results;
        for (double angle : angles)
        {
            cds::RotationMatrix matrix = rotationTo(dihedral, constants::toRadians(angle));
            applyMatrix(input.graph, input.bounds, movingAtoms, matrix, bounds);
            cds::Overlap overlaps = cds::overlapVectorSum(cds::CountOverlappingAtoms(
                potential, overlapProperties, input.graph, bounds, input.residueWeights, input.atomElements,
                input.atomIncluded, input.bonds, fixedResidueIndices, movingResidueIndices));

            cds::AngleOverlap current {
                overlaps, cds::AngleWithMetadata {angle, anglePreference, metadataIndex}
            };
            if (overlaps.count <= 0.0)
            {
                // requires that the angles are sorted in order of closest to preference
                return current;
            }
            else
            {
                results.push_back(current);
            }
        }
        size_t index = bestOverlapResultIndex(results);
        return results[index];
    }
} // namespace

size_t cds::bestOverlapResultIndex(const std::vector<AngleOverlap>& results)
{
    auto differenceFromPreference = [](AngleOverlap& a)
    {
        return std::abs(a.angle.value - a.angle.preference);
    };
    size_t bestIndex = 0;
    for (size_t n = 1; n < results.size(); n++)
    {
        auto a            = results[n];
        auto b            = results[bestIndex];
        int comp          = compareOverlaps(a.overlaps, b.overlaps);
        bool sameMetadata = (a.angle.metadataIndex == b.angle.metadataIndex);
        if ((comp < 0) || (sameMetadata && (comp == 0) && differenceFromPreference(a) < differenceFromPreference(b)))
        {
            bestIndex = n;
        }
    }
    return bestIndex;
}

cds::DihedralRotationDataContainer
cds::dihedralRotationInputData(double overlapTolerance, RotatableDihedral& dihedral, const GraphIndexData& indices,
                               const std::array<ResiduesWithOverlapWeight, 2>& residueSets)
{
    assembly::Graph graph  = createAssemblyGraph(indices, std::vector<bool>(indices.atoms.size(), true));
    auto movingAtomSpheres = atomCoordinatesWithRadii(dihedral.movingAtoms);
    auto movingAtomBounds  = boundingSphere(movingAtomSpheres);
    Coordinate pointA      = dihedral.atoms[1]->coordinate();
    Coordinate pointB      = dihedral.atoms[2]->coordinate();
    auto movementBounds    = boundingSphereCenteredOnLine(movingAtomBounds, pointA, pointB);

    const std::vector<cds::Residue*>& residues = indices.residues;
    std::vector<bool> atomMoving               = toAtomMoving(indices.atoms, dihedral.movingAtoms);
    std::vector<bool> atomIncluded(graph.atomCount, true);
    std::vector<bool> residueMoving(graph.residueCount, false);
    for (size_t n = 0; n < graph.atomResidue.size(); n++)
    {
        if (atomMoving[n])
        {
            residueMoving[graph.atomResidue[n]] = true;
        }
    }
    std::vector<size_t> residueSetA = toResidueIndices(residues, residueSets[0].residues);
    std::vector<size_t> residueSetB = toResidueIndices(residues, residueSets[1].residues);
    std::vector<double> residueWeights(graph.residueCount, 0.0);
    for (size_t n = 0; n < residueSetA.size(); n++)
    {
        residueWeights[residueSetA[n]] = residueSets[0].weights[n];
    }
    for (size_t n = 0; n < residueSetB.size(); n++)
    {
        residueWeights[residueSetB[n]] = residueSets[1].weights[n];
    }
    std::array<std::vector<size_t>, 2> residueIndices =
        dihedral.isBranchingLinkage ? branchedResidueSets(residueMoving, residueSetA, residueSetB)
                                    : std::array<std::vector<size_t>, 2> {residueSetA, residueSetB};

    std::vector<Atom*> atomsA = codeUtils::indicesToValues(indices.atoms, residueAtoms(graph, residueIndices[0][0]));
    std::vector<Atom*> atomsB = codeUtils::indicesToValues(indices.atoms, residueAtoms(graph, residueIndices[1][0]));
    auto residueBond          = bondedAtomPair(atomsA, atomsB);
    auto bonds                = std::vector<cds::BondedResidueOverlapInput> {
        {{residueIndices[0][0], residueIndices[1][0]},
         {atomsBondedTo(residueBond[0], atomsA), atomsBondedTo(residueBond[1], atomsB)}}
    };
    std::vector<Sphere> atomBounds = atomCoordinatesWithRadii(indices.atoms);
    std::vector<Sphere> residueBounds;
    residueBounds.reserve(graph.residueCount);
    for (size_t n = 0; n < graph.residueCount; n++)
    {
        residueBounds.push_back(boundingSphere(codeUtils::indicesToValues(atomBounds, residueAtoms(graph, n))));
    }
    std::vector<Sphere> moleculeBounds;
    moleculeBounds.reserve(graph.moleculeCount);
    for (size_t n = 0; n < graph.moleculeCount; n++)
    {
        moleculeBounds.push_back(boundingSphere(codeUtils::indicesToValues(residueBounds, moleculeResidues(graph, n))));
    }
    assembly::Bounds bounds {atomBounds, residueBounds, moleculeBounds};
    return {
        graph,
        bounds,
        atomMoving,
        atomIncluded,
        cds::atomElementEnums(indices.atoms),
        residueWeights,
        {residueIndices[0],
          cds::intersectingIndices(overlapTolerance, movementBounds, residueBounds, residueIndices[1])},
        bonds
    };
}

cds::OverlapState cds::wiggleUsingRotamers(const MolecularMetadata::PotentialTable& potential,
                                           cds::OverlapProperties overlapProperties, SearchAngles searchAngles,
                                           const cds::DihedralCoordinates coordinates,
                                           const std::vector<size_t>& indices, const DihedralAngleDataVector& rotamers,
                                           const AngleSearchPreference& preference,
                                           const cds::DihedralRotationData& input)
{
    // residue sets can be empty if overlap tolerance is high
    // no overlaps are possible in this case, so we can
    // return preferred angle immediately to avoid segfault later
    if (input.residueIndices[0].empty() || input.residueIndices[1].empty())
    {
        size_t index = preference.metadataOrder[0];
        double angle = preference.angles[index];
        return {
            cds::Overlap {0.0, 0.0},
             {angle, angle, index},
             input.bounds
        };
    }
    else
    {
        auto resultState = [&](const AngleOverlap best)
        {
            assembly::Bounds bounds    = input.bounds;
            cds::RotationMatrix matrix = rotationTo(coordinates, constants::toRadians(best.angle.value));
            applyMatrix(input.graph, input.bounds, codeUtils::boolsToIndices(input.atomMoving), matrix, bounds);
            return OverlapState {best.overlaps, best.angle, bounds};
        };
        std::vector<AngleOverlap> results;
        for (size_t n : preference.metadataOrder)
        {
            double angle     = preference.angles[n];
            double deviation = preference.deviation;
            AngleOverlap best =
                WiggleAnglesOverlaps(potential, overlapProperties, coordinates, indices[n], preference.angles[n],
                                     searchAngles(rotamers[n], angle, deviation), input);
            // found something with no overlaps
            // if metadata and angles are sorted in order of preference, we can quit here
            if (best.overlaps.count <= 0.0)
            {
                return resultState(best);
            }
            results.push_back(best);
        }
        size_t index = bestOverlapResultIndex(results);
        return resultState(results[index]);
    }
}

void cds::simpleWiggleCurrentRotamers(const MolecularMetadata::PotentialTable& potential,
                                      cds::OverlapProperties overlapProperties, SearchAngles searchAngles,
                                      std::vector<RotatableDihedral>& dihedrals,
                                      const std::vector<DihedralAngleDataVector>& metadata,
                                      const std::vector<AngleSearchPreference>& preference,
                                      const GraphIndexData& indices,
                                      const std::array<ResiduesWithOverlapWeight, 2>& residues)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        cds::RotatableDihedral& dihedral = dihedrals[n];
        std::vector<size_t> index {dihedral.currentMetadataIndex};
        std::array<Coordinate, 4> coordinates = dihedralCoordinates(dihedral);
        DihedralRotationDataContainer input =
            dihedralRotationInputData(overlapProperties.tolerance, dihedral, indices, residues);
        DihedralRotationData inputPointers {
            input.graph,
            input.bounds,
            input.atomMoving,
            input.atomIncluded,
            input.atomElements,
            input.residueWeights,
            {input.residueIndices[0], input.residueIndices[1]},
            input.bonds
        };
        OverlapState best = wiggleUsingRotamers(potential, overlapProperties, searchAngles, coordinates, index,
                                                metadata[n], preference[n], inputPointers);
        for (size_t n = 0; n < input.graph.atomCount; n++)
        {
            if (input.atomMoving[n])
            {
                indices.atoms[n]->setCoordinate(best.bounds.atoms[n].center);
            }
        }
    }
}

std::vector<double> cds::evenlySpacedAngles(double preference, double lowerDeviation, double upperDeviation,
                                            double increment)
{
    auto closerToPreference = [&preference](double a, double b)
    {
        return std::abs(a - preference) < std::abs(b - preference);
    };
    auto lowerRange            = evenlySpaced(preference - lowerDeviation, preference, increment);
    auto upperRange            = evenlySpaced(preference + upperDeviation, preference, increment);
    std::vector<double> result = codeUtils::vectorAppend(lowerRange, upperRange);
    result.push_back(preference);
    // sorted angles enable early return from angle search when 0 overlaps found
    std::sort(result.begin(), result.end(), closerToPreference);
    return result;
}

std::vector<cds::AngleSearchPreference> cds::angleSearchPreference(double deviation,
                                                                   const ResidueLinkageShapePreference& preference)
{
    std::function<std::vector<cds::AngleSearchPreference>(const ConformerShapePreference&)> onConformer =
        [&](const ConformerShapePreference& pref)
    {
        if (pref.metadataOrder.size() != 1)
        {
            throw std::runtime_error(
                "cannot construct search preference for conformer with multiple available rotamers");
        }
        std::vector<AngleSearchPreference> result;
        result.reserve(pref.angles.size());
        for (auto& a : pref.angles)
        {
            result.push_back({deviation, a, pref.metadataOrder});
        }
        return result;
    };
    std::function<std::vector<cds::AngleSearchPreference>(const PermutationShapePreference&)> onPermutation =
        [&](const PermutationShapePreference& pref)
    {
        std::vector<AngleSearchPreference> result;
        result.reserve(pref.angles.size());
        for (size_t n = 0; n < pref.angles.size(); n++)
        {
            result.push_back({deviation, pref.angles[n], pref.metadataOrder[n]});
        }
        return result;
    };
    return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
}

std::vector<std::vector<cds::AngleSearchPreference>>
cds::angleSearchPreference(double deviation, const std::vector<ResidueLinkageShapePreference>& preferences)
{
    std::vector<std::vector<AngleSearchPreference>> result;
    result.reserve(preferences.size());
    for (auto& a : preferences)
    {
        result.push_back(angleSearchPreference(deviation, a));
    }
    return result;
}
