#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"

#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
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

    std::vector<cds::Atom*> movingAtomsWithinSet(cds::Atom* dihedralFixed, cds::Atom* dihedralMoving,
                                                 const std::vector<cds::Atom*> set)
    {
        std::vector<cds::Atom*> result;
        std::vector<cds::Atom*> queued;
        auto shouldAdd = [&](cds::Atom* atom)
        {
            return codeUtils::contains(set, atom) && !codeUtils::contains(result, atom) &&
                   !codeUtils::contains(queued, atom);
        };
        result.push_back(dihedralFixed);
        queued.push_back(dihedralMoving);
        while (!queued.empty())
        {
            auto next = queued.back();
            result.push_back(next);
            queued.pop_back();
            for (auto& a : next->getNeighbors())
            {
                if (shouldAdd(a))
                {
                    queued.push_back(a);
                }
            }
        }
        result.erase(result.begin());
        return result;
    }

    std::vector<bool> movingAtomsWithinResidue(const std::vector<cds::Atom*> movingAtoms,
                                               const std::vector<cds::Atom*> residueAtoms)
    {
        std::vector<bool> result;
        for (auto& a : residueAtoms)
        {
            result.push_back(codeUtils::contains(movingAtoms, a));
        }
        return result;
    }

    std::array<cds::ResiduesWithOverlapWeight, 2>
    branchedResidueSets(const std::vector<cds::CoordinateReference>& movingCoordinates,
                        const std::array<cds::ResiduesWithOverlapWeight, 2>& input)
    {
        auto residueContainsMovingAtom = [&](cds::Residue* res)
        {
            return codeUtils::contains(movingCoordinates, res->getAtoms()[0]->coordinateReference());
        };
        auto& setA                         = input[0];
        auto& setB                         = input[1];
        std::vector<cds::Residue*> newSetA = setB.residues;
        std::vector<cds::Residue*> newSetB {setA.residues[0]};
        std::vector<double> newWeightsA = input[1].weights;
        std::vector<double> newWeightsB {input[0].weights[0]};
        for (size_t n = 1; n < setA.residues.size(); n++)
        {
            auto res      = setA.residues[n];
            double weight = setA.weights[n];
            if (residueContainsMovingAtom(res))
            {
                newSetB.push_back(res);
                newWeightsB.push_back(weight);
            }
            else
            {
                newSetA.push_back(res);
                newWeightsA.push_back(weight);
            }
        }
        return {
            cds::ResiduesWithOverlapWeight {newSetA, newWeightsA},
            cds::ResiduesWithOverlapWeight {newSetB, newWeightsB}
        };
    }

    cds::DihedralRotationData toRotationData(const cds::Sphere bounds, const std::vector<cds::Atom*>& movingAtoms,
                                             const std::vector<cds::Residue*>& residues,
                                             const std::vector<double>& residueWeights,
                                             const std::vector<bool>& firstResidueBondedAtoms)
    {
        size_t atomCount = 0;
        for (auto res : residues)
        {
            atomCount += res->atomCount();
        }
        std::vector<cds::Sphere> coordinates;
        coordinates.reserve(atomCount);
        std::vector<std::vector<size_t>> residueAtoms;
        residueAtoms.reserve(residues.size());
        size_t currentAtom = 0;
        for (auto& res : residues)
        {
            auto& atomsRef = res->getAtomsReference();
            residueAtoms.push_back(codeUtils::indexVectorWithOffset(currentAtom, atomsRef));
            for (auto& atomPtr : atomsRef)
            {
                coordinates.push_back(coordinateWithRadius(atomPtr.get()));
            }
            currentAtom += atomsRef.size();
        }
        std::vector<std::vector<size_t>> withinRangeResidueAtoms;
        withinRangeResidueAtoms.reserve(residues.size());
        std::vector<cds::Sphere> withinRangeSpheres;
        withinRangeSpheres.reserve(residues.size());
        std::vector<double> withinRangeWeights;
        withinRangeWeights.reserve(residues.size());
        std::vector<cds::Sphere> residuePoints;
        for (size_t n = 0; n < residueAtoms.size(); n++)
        {
            // always include first residue, we assume that it's bonded with first residue of the other set
            bool within                  = n == 0;
            std::vector<size_t>& indices = residueAtoms[n];
            for (size_t index : indices)
            {
                within = within || cds::spheresOverlap(constants::overlapTolerance, bounds, coordinates[index]);
                if (within)
                {
                    break;
                }
            }
            if (within)
            {
                residuePoints.clear();
                for (size_t index : indices)
                {
                    residuePoints.push_back(coordinates[index]);
                }
                withinRangeSpheres.push_back(boundingSphere(residuePoints));
                withinRangeWeights.push_back(residueWeights[n]);
                withinRangeResidueAtoms.push_back(indices);
            }
        }

        const auto& firstAtoms   = residues[0]->getAtoms();
        std::vector<bool> moving = movingAtomsWithinResidue(movingAtoms, firstAtoms);
        return {coordinates, withinRangeSpheres,     withinRangeResidueAtoms, withinRangeWeights,
                moving,      firstResidueBondedAtoms};
    }

    void moveFirstResidueCoords(size_t atomOffset, size_t residueOffset, const cds::RotationMatrix matrix,
                                const cds::DihedralRotationData& input, std::vector<cds::Sphere>& coordinates,
                                std::vector<cds::Sphere>& spheres)
    {
        std::vector<cds::Sphere> firstCoords;
        const std::vector<size_t>& indices = input.residueAtoms[0];
        firstCoords.reserve(indices.size());
        for (size_t index : indices)
        {
            const Coordinate& center = input.coordinates[index].center;
            coordinates[index + atomOffset].center =
                input.firstResidueCoordinateMoving[index] ? matrix * center : center;
            firstCoords.push_back(coordinates[index + atomOffset]);
        }
        spheres[residueOffset] = cds::boundingSphere(firstCoords);
    }

    cds::AngleOverlap WiggleAnglesOverlaps(const cds::DihedralCoordinates dihedral, size_t metadataIndex,
                                           double anglePreference, std::vector<double> angles,
                                           const std::array<cds::DihedralRotationData, 2>& input)
    {
        auto& fixedInput                           = input[0];
        auto& movingInput                          = input[1];
        // copies of input vectors used for updates during looping
        std::vector<cds::Sphere> fixedCoordinates  = fixedInput.coordinates;
        std::vector<cds::Sphere> fixedSpheres      = fixedInput.boundingSpheres;
        std::vector<cds::Sphere> movingCoordinates = movingInput.coordinates;
        std::vector<cds::Sphere> movingSpheres     = movingInput.boundingSpheres;
        size_t atomOffset                          = fixedCoordinates.size();
        size_t residueOffset                       = fixedSpheres.size();
        std::vector<cds::Sphere> atomBounds        = codeUtils::vectorAppend(fixedCoordinates, movingCoordinates);
        std::vector<cds::Sphere> residueBounds     = codeUtils::vectorAppend(fixedSpheres, movingSpheres);
        std::vector<double> weights = codeUtils::vectorAppend(fixedInput.residueWeights, movingInput.residueWeights);
        std::vector<size_t> fixedResidueIndices = codeUtils::indexVector(fixedInput.residueAtoms);
        std::vector<size_t> movingResidueIndices =
            codeUtils::indexVectorWithOffset(residueOffset, movingInput.residueAtoms);

        auto movingAtomIndices = movingInput.residueAtoms;
        for (size_t n = 0; n < movingAtomIndices.size(); n++)
        {
            movingAtomIndices[n] = codeUtils::offsetIndices(atomOffset, movingAtomIndices[n]);
        }
        std::vector<std::vector<size_t>> residueAtoms =
            codeUtils::vectorAppend(fixedInput.residueAtoms, movingAtomIndices);
        std::vector<cds::BondedResidueOverlapInput> bonds {
            {{fixedResidueIndices[0], movingResidueIndices[0]},
             {fixedInput.firstResidueBondedAtoms, movingInput.firstResidueBondedAtoms}}
        };
        std::vector<cds::AngleOverlap> results;
        for (double angle : angles)
        {
            auto matrix = rotationTo(dihedral, constants::toRadians(angle));
            moveFirstResidueCoords(0, 0, matrix, fixedInput, atomBounds, residueBounds);
            moveFirstResidueCoords(atomOffset, residueOffset, matrix, movingInput, atomBounds, residueBounds);
            for (size_t n = 1; n < movingInput.residueAtoms.size(); n++)
            {
                for (size_t index : movingInput.residueAtoms[n])
                {
                    atomBounds[index + atomOffset].center = matrix * movingInput.coordinates[index].center;
                }
            }
            for (size_t n = 1; n < movingResidueIndices.size(); n++)
            {
                residueBounds[movingResidueIndices[n]].center = matrix * movingInput.boundingSpheres[n].center;
            }
            cds::Overlap overlaps = cds::CountOverlappingAtoms({atomBounds, residueBounds, residueAtoms, weights},
                                                               bonds, fixedResidueIndices, movingResidueIndices);

            results.push_back({
                overlaps, cds::AngleWithMetadata {angle, anglePreference, metadataIndex}
            });
        }
        return bestOverlapResult(results);
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

cds::AngleOverlap cds::bestOverlapResult(const std::vector<AngleOverlap>& results)
{
    return results[bestOverlapResultIndex(results)];
}

std::array<cds::DihedralRotationData, 2>
cds::dihedralRotationInputData(RotatableDihedral& dihedral, const std::array<ResiduesWithOverlapWeight, 2>& residues)
{
    auto& atoms                      = dihedral.atoms;
    auto movingCoordinates           = atomCoordinateReferences(dihedral.movingAtoms);
    auto dihedralResiduesMovingAtoms = movingAtomsWithinSet(
        atoms[2], atoms[1],
        codeUtils::vectorAppend(residues[0].residues[0]->getAtoms(), residues[1].residues[0]->getAtoms()));

    auto movingAtomSpheres  = atomCoordinatesWithRadii(dihedral.movingAtoms);
    auto movingAtomBounds   = boundingSphere(movingAtomSpheres);
    Coordinate origin       = atoms[1]->coordinate();
    Coordinate axis         = atoms[2]->coordinate() - origin;
    auto closestPointOnAxis = origin + projection(movingAtomBounds.center - origin, axis);
    double distanceToAxis   = length(closestPointOnAxis - movingAtomBounds.center);
    auto movementBounds     = Sphere {movingAtomBounds.radius + distanceToAxis, closestPointOnAxis};

    auto rotationData = [&](const cds::ResiduesWithOverlapWeight& set, const std::vector<bool>& firstResidueBondedAtoms)
    {
        return toRotationData(movementBounds, dihedralResiduesMovingAtoms, set.residues, set.weights,
                              firstResidueBondedAtoms);
    };

    auto inputSets   = dihedral.isBranchingLinkage ? branchedResidueSets(movingCoordinates, residues) : residues;
    auto atomsA      = inputSets[0].residues[0]->getAtoms();
    auto atomsB      = inputSets[1].residues[0]->getAtoms();
    auto residueBond = bondedAtomPair(atomsA, atomsB);
    return {rotationData(inputSets[0], atomsBondedTo(residueBond[0], atomsA)),
            rotationData(inputSets[1], atomsBondedTo(residueBond[1], atomsB))};
}

cds::AngleOverlap cds::wiggleUsingRotamers(SearchAngles searchAngles, const cds::DihedralCoordinates coordinates,
                                           const std::vector<size_t>& indices, const DihedralAngleDataVector& rotamers,
                                           const AngleSearchPreference& preference,
                                           const std::array<cds::DihedralRotationData, 2>& input)
{
    std::vector<AngleOverlap> results;
    for (size_t n : preference.metadataOrder)
    {
        double angle      = preference.angles[n];
        double deviation  = preference.deviation;
        AngleOverlap best = WiggleAnglesOverlaps(coordinates, indices[n], preference.angles[n],
                                                 searchAngles(rotamers[n], angle, deviation), input);
        // found something with no overlaps
        // if metadata and angles are sorted in order of preference, we can quit here
        if (best.overlaps.count <= 0.0)
        {
            return best;
        }
        results.push_back(best);
    }

    return bestOverlapResult(results);
}

void cds::simpleWiggleCurrentRotamers(SearchAngles searchAngles, std::vector<RotatableDihedral>& dihedrals,
                                      const std::vector<DihedralAngleDataVector>& metadata,
                                      const std::vector<AngleSearchPreference>& preference,
                                      const std::array<ResiduesWithOverlapWeight, 2>& residues)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        auto& dihedral = dihedrals[n];
        std::vector<size_t> index {dihedral.currentMetadataIndex};
        auto coordinates = dihedralCoordinates(dihedral);
        auto input       = dihedralRotationInputData(dihedral, residues);
        auto best        = wiggleUsingRotamers(searchAngles, coordinates, index, metadata[n], preference[n], input);
        setDihedralAngle(dihedral, best.angle);
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
    std::sort(result.begin(), result.end(), closerToPreference);
    return result;
}

std::vector<cds::AngleSearchPreference> cds::angleSearchPreference(double deviation,
                                                                   const ResidueLinkageShapePreference& preference)
{
    if (std::holds_alternative<ConformerShapePreference>(preference))
    {
        auto pref = std::get<ConformerShapePreference>(preference);
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
    }
    else if (std::holds_alternative<PermutationShapePreference>(preference))
    {

        auto pref = std::get<PermutationShapePreference>(preference);
        std::vector<AngleSearchPreference> result;
        result.reserve(pref.angles.size());
        for (size_t n = 0; n < pref.angles.size(); n++)
        {
            result.push_back({deviation, pref.angles[n], pref.metadataOrder[n]});
        }
        return result;
    }
    else
    {
        throw std::runtime_error("unhandled linkage shape preference in cds::angleSearchPreference");
    }
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
