#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"

#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <numeric>
#include <vector>

using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;

namespace
{
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
    branchedResidueSets(const std::vector<Coordinate*>& movingCoordinates,
                        const std::array<cds::ResiduesWithOverlapWeight, 2>& input)
    {
        auto residueContainsMovingAtom = [&](cds::Residue* res)
        {
            return codeUtils::contains(movingCoordinates, res->getAtoms()[0]->getCoordinate());
        };
        auto& setA                         = input[0];
        auto& setB                         = input[1];
        std::vector<cds::Residue*> newSetA = setB.residues;
        std::vector<cds::Residue*> newSetB {setA.residues[0]};
        std::vector<double> newWeightsA = input[1].weights;
        std::vector<double> newWeightsB {input[0].weights[0]};
        for (size_t n = 1; n < setA.residues.size(); n++)
        {
            auto res    = setA.residues[n];
            uint weight = setA.weights[n];
            if (residueContainsMovingAtom(res))
            {
                newSetB.push_back(res);
                newWeightsB.push_back(weight);
            }
            else
            {
                newSetA.push_back(res);
                newWeightsB.push_back(weight);
            }
        }
        return {
            cds::ResiduesWithOverlapWeight {newSetA, newWeightsA},
            cds::ResiduesWithOverlapWeight {newSetB, newWeightsB}
        };
    }

    cds::DihedralRotationData toRotationData(const cds::Sphere bounds, const std::vector<cds::Atom*>& movingAtoms,
                                             const std::vector<cds::Residue*>& residues,
                                             const std::vector<double>& residueWeights)
    {
        size_t atomCount = 0;
        for (auto res : residues)
        {
            atomCount += res->atomCount();
        }
        std::vector<cds::Sphere> coordinates;
        coordinates.reserve(atomCount);
        std::vector<std::pair<size_t, size_t>> residueAtoms;
        residueAtoms.reserve(residues.size());
        size_t currentAtom = 0;
        for (auto& res : residues)
        {
            size_t startAtom = currentAtom;
            for (const auto& atomPtr : res->getAtomsReference())
            {
                coordinates.push_back(coordinateWithRadius(atomPtr.get()));
            }
            currentAtom += res->atomCount();
            residueAtoms.push_back({startAtom, currentAtom});
        }
        std::vector<std::pair<size_t, size_t>> withinRangeResidueAtoms;
        withinRangeResidueAtoms.reserve(residues.size());
        std::vector<cds::Sphere> withinRangeSpheres;
        withinRangeSpheres.reserve(residues.size());
        std::vector<double> withinRangeWeights;
        withinRangeWeights.reserve(residues.size());
        std::vector<cds::Sphere> residuePoints;
        for (size_t n = 0; n < residueAtoms.size(); n++)
        {
            auto range  = residueAtoms[n];
            bool within = false;
            for (size_t k = range.first; k < range.second; k++)
            {
                within = cds::spheresOverlap(constants::overlapTolerance, bounds, coordinates[k]);
                if (within)
                {
                    break;
                }
            }
            if (within)
            {
                residuePoints.clear();
                residuePoints.insert(residuePoints.end(), coordinates.begin() + range.first,
                                     coordinates.begin() + range.second);
                withinRangeSpheres.push_back(boundingSphere(residuePoints));
                withinRangeWeights.push_back(residueWeights[n]);
                withinRangeResidueAtoms.push_back(range);
            }
        }

        const auto& firstAtoms   = residues[0]->getAtoms();
        std::vector<bool> moving = movingAtomsWithinResidue(movingAtoms, firstAtoms);
        return {coordinates, withinRangeSpheres, withinRangeResidueAtoms, withinRangeWeights, moving};
    }

    double maxDistanceFrom(const Coordinate& pt, const std::vector<cds::Coordinate*>& coords)
    {
        double maxSquare = 0.0;
        for (auto& a : coords)
        {
            maxSquare = std::max(maxSquare, squaredDistance(pt, *a));
        }
        return std::sqrt(maxSquare);
    }

    std::vector<double> wiggleAngles(cds::Bounds bounds, double approximateIncrement)
    {
        double range     = bounds.upper - bounds.lower;
        int steps        = std::ceil(range / approximateIncrement);
        double increment = range / steps;
        std::vector<double> result;
        result.reserve(steps + 1);
        for (int k = 0; k < steps + 1; k++)
        {
            result.push_back(bounds.lower + k * increment);
        }
        return result;
    }

    void moveFirstResidueCoords(const cds::RotationMatrix matrix, const cds::DihedralRotationData& input,
                                std::vector<cds::Sphere>& coordinates, std::vector<cds::Sphere>& spheres)
    {
        std::vector<cds::Sphere> firstCoords;
        auto range = input.residueAtoms[0];
        firstCoords.reserve(range.second - range.first);
        for (size_t n = range.first; n < range.second; n++)
        {
            auto& center          = input.coordinates[n].center;
            coordinates[n].center = input.firstResidueCoordinateMoving[n] ? matrix * center : center;
            firstCoords.push_back(coordinates[n]);
        }
        spheres[0] = cds::boundingSphere(firstCoords);
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
        std::vector<cds::AngleOverlap> results;
        for (double angle : angles)
        {
            auto matrix = rotationTo(dihedral, constants::toRadians(angle));
            moveFirstResidueCoords(matrix, fixedInput, fixedCoordinates, fixedSpheres);
            moveFirstResidueCoords(matrix, movingInput, movingCoordinates, movingSpheres);
            for (size_t n = movingInput.residueAtoms[0].second; n < movingCoordinates.size(); n++)
            {
                movingCoordinates[n].center = matrix * movingInput.coordinates[n].center;
            }
            for (size_t n = 1; n < movingSpheres.size(); n++)
            {
                movingSpheres[n].center = matrix * movingInput.boundingSpheres[n].center;
            }
            cds::Overlap overlaps = cds::CountOverlappingAtoms(
                {fixedCoordinates, fixedSpheres, fixedInput.residueAtoms, fixedInput.residueWeights},
                {movingCoordinates, movingSpheres, movingInput.residueAtoms, movingInput.residueWeights});

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
        auto a   = results[n];
        auto b   = results[bestIndex];
        int comp = compareOverlaps(a.overlaps, b.overlaps);
        if ((comp < 0) || ((comp == 0) && differenceFromPreference(a) < differenceFromPreference(b)))
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
    auto& movingCoordinates          = dihedral.movingCoordinates;
    auto dihedralResiduesMovingAtoms = movingAtomsWithinSet(
        atoms[2], atoms[1],
        codeUtils::vectorAppend(residues[0].residues[0]->getAtoms(), residues[1].residues[0]->getAtoms()));

    auto centerPoint      = *atoms[1]->getCoordinate();
    double maxDistance    = maxDistanceFrom(centerPoint, movingCoordinates);
    double cutoffDistance = maxDistance + constants::maxCutOff;

    auto rotationData = [&](const cds::ResiduesWithOverlapWeight& set)
    {
        return toRotationData(Sphere {cutoffDistance, centerPoint}, dihedralResiduesMovingAtoms, set.residues,
                              set.weights);
    };

    auto inputSets = dihedral.isBranchingLinkage ? branchedResidueSets(movingCoordinates, residues) : residues;
    return {rotationData(inputSets[0]), rotationData(inputSets[1])};
}

cds::AngleOverlap cds::wiggleWithinRangesDistanceCheck(RotatableDihedral& dihedral,
                                                       const AngleSearchPreference& preference,
                                                       std::vector<cds::Atom*>& overlapAtomSet1,
                                                       std::vector<cds::Atom*>& overlapAtomSet2, size_t metadataIndex,
                                                       std::vector<double> angles)
{
    size_t currentMetadataIndex = dihedral.currentMetadataIndex;
    AngleWithMetadata initial   = {cds::angle(dihedralCoordinates(dihedral)), preference.angles[currentMetadataIndex],
                                   currentMetadataIndex};
    std::vector<cds::AngleOverlap> results;
    for (double angle : angles)
    {
        AngleWithMetadata angleWithMetadata = {angle, preference.angles[metadataIndex], metadataIndex};
        setDihedralAngle(dihedral, angleWithMetadata);
        cds::Overlap overlaps = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);

        results.push_back({overlaps, angleWithMetadata});
    }
    setDihedralAngle(dihedral, initial);
    return bestOverlapResult(results);
}

// User requested gg, this prevents flipping into gt like the above would do. i.e. cb won't want a flip, gp would.
cds::AngleOverlap cds::wiggleWithinCurrentRotamer(cds::RotatableDihedral& dihedral,
                                                  const DihedralAngleDataVector& metadataVector,
                                                  const AngleSearchPreference& preference,
                                                  std::vector<cds::Atom*>& overlapAtomSet1,
                                                  std::vector<cds::Atom*>& overlapAtomSet2, double angleIncrement)
{
    size_t currentMetadataIndex = dihedral.currentMetadataIndex;
    if (preference.metadataOrder.size() != 1 || preference.metadataOrder[0] != currentMetadataIndex)
    {
        throw std::runtime_error("wrong preference for current rotamer");
    }
    Bounds bounds = angleBounds(metadataVector[currentMetadataIndex]);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    return wiggleWithinRangesDistanceCheck(dihedral, preference, overlapAtomSet1, overlapAtomSet2, currentMetadataIndex,
                                           wiggleAngles(bounds, angleIncrement));
}

cds::AngleOverlap cds::wiggleUsingRotamers(const cds::DihedralCoordinates coordinates,
                                           const std::vector<size_t>& indices, const DihedralAngleDataVector& rotamers,
                                           const AngleSearchPreference& preference, double angleIncrement,
                                           const std::array<cds::DihedralRotationData, 2>& input)
{
    std::vector<AngleOverlap> results;
    for (size_t n : preference.metadataOrder)
    {
        Bounds bounds     = angleBounds(rotamers[n]);
        AngleOverlap best = WiggleAnglesOverlaps(coordinates, indices[n], preference.angles[n],
                                                 wiggleAngles(bounds, angleIncrement), input);
        results.push_back(best);
    }

    return bestOverlapResult(results);
}

void cds::simpleWiggleCurrentRotamers(std::vector<RotatableDihedral>& dihedrals,
                                      const std::vector<DihedralAngleDataVector>& metadata,
                                      const std::vector<AngleSearchPreference>& preferences,
                                      std::vector<cds::Atom*>& overlapAtomSet1,
                                      std::vector<cds::Atom*>& overlapAtomSet2, double angleIncrement)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        auto& dihedral = dihedrals[n];
        auto best = wiggleWithinCurrentRotamer(dihedral, metadata[n], preferences[n], overlapAtomSet1, overlapAtomSet2,
                                               angleIncrement);
        setDihedralAngle(dihedral, best.angle);
    }
}

void cds::simpleWiggleCurrentRotamers(std::vector<RotatableDihedral>& dihedrals,
                                      const std::vector<DihedralAngleDataVector>& metadata,
                                      const std::vector<AngleSearchPreference>& preference,
                                      const std::array<ResiduesWithOverlapWeight, 2>& residues, double angleIncrement)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        auto& dihedral = dihedrals[n];
        std::vector<size_t> index {dihedral.currentMetadataIndex};
        auto coordinates = dihedralCoordinates(dihedral);
        auto input       = dihedralRotationInputData(dihedral, residues);
        auto best        = wiggleUsingRotamers(coordinates, index, metadata[n], preference[n], angleIncrement, input);
        setDihedralAngle(dihedral, best.angle);
    }
}

std::vector<cds::AngleSearchPreference> cds::angleSearchPreference(const ResidueLinkageShapePreference& preference)
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
            result.push_back({a, pref.metadataOrder});
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
            result.push_back({pref.angles[n], pref.metadataOrder[n]});
        }
        return result;
    }
    else
    {
        throw std::runtime_error("unhandled linkage shape preference in cds::angleSearchPreference");
    }
}

std::vector<std::vector<cds::AngleSearchPreference>>
cds::angleSearchPreference(const std::vector<ResidueLinkageShapePreference>& preferences)
{
    std::vector<std::vector<AngleSearchPreference>> result;
    result.reserve(preferences.size());
    for (auto& a : preferences)
    {
        result.push_back(angleSearchPreference(a));
    }
    return result;
}
