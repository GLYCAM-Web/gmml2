#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"

#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"
#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <algorithm>
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

    std::array<std::vector<cds::Residue*>, 2>
    branchedResidueSets(const std::vector<Coordinate*>& movingCoordinates,
                        const std::array<std::vector<cds::Residue*>, 2>& residues)
    {
        auto residueContainsMovingAtom = [&](cds::Residue* res)
        {
            return codeUtils::contains(movingCoordinates, res->getAtoms()[0]->getCoordinate());
        };
        auto& setA                         = residues[0];
        auto& setB                         = residues[1];
        std::vector<cds::Residue*> newSetA = setB;
        std::vector<cds::Residue*> newSetB {setA[0]};
        for (size_t n = 1; n < setA.size(); n++)
        {
            auto res = setA[n];
            if (residueContainsMovingAtom(res))
            {
                newSetB.push_back(res);
            }
            else
            {
                newSetA.push_back(res);
            }
        }
        return {newSetA, newSetB};
    }

    cds::DihedralRotationData toRotationData(const cds::Sphere bounds, const std::vector<cds::Atom*> movingAtoms,
                                             const std::vector<cds::Residue*> residues)
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
                withinRangeResidueAtoms.push_back(range);
            }
        }

        const auto& firstAtoms   = residues[0]->getAtoms();
        std::vector<bool> moving = movingAtomsWithinResidue(movingAtoms, firstAtoms);
        return {coordinates, withinRangeSpheres, withinRangeResidueAtoms, moving};
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

    std::vector<double> wiggleAngles(cds::Bounds bounds, int approximateIncrement)
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

    cds::AngleOverlap bestOverlapResult(const std::vector<cds::AngleOverlap>& results)
    {
        auto differenceFromDefault = [](cds::AngleOverlap& a)
        {
            return std::abs(a.angle.value - a.angle.metadata.default_angle_value_);
        };
        int best = 0;
        for (size_t n = 1; n < results.size(); n++)
        {
            auto a   = results[n];
            auto b   = results[best];
            int comp = cds::compareOverlaps(a.overlaps, b.overlaps);
            if ((comp < 0) || ((comp == 0) && differenceFromDefault(a) < differenceFromDefault(b)))
            {
                best = n;
            }
        }
        return results[best];
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

    cds::AngleOverlap WiggleAnglesOverlaps(const cds::DihedralCoordinates dihedral, const DihedralAngleData& metadata,
                                           std::vector<double> angles,
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
            auto matrix = rotationTo(dihedral, constants::degree2Radian(angle));
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
            cds::Overlap overlaps =
                cds::CountOverlappingAtoms({fixedCoordinates, fixedSpheres, fixedInput.residueAtoms},
                                           {movingCoordinates, movingSpheres, movingInput.residueAtoms});

            results.push_back({
                overlaps, cds::AngleWithMetadata {angle, metadata}
            });
        }
        return bestOverlapResult(results);
    }
} // namespace

std::array<cds::DihedralRotationData, 2>
cds::dihedralRotationInputData(RotatableDihedral& dihedral, const std::array<std::vector<Residue*>, 2>& residues)
{
    auto& atoms                      = dihedral.atoms;
    auto& movingCoordinates          = dihedral.movingCoordinates;
    auto dihedralResiduesMovingAtoms = movingAtomsWithinSet(
        atoms[2], atoms[1], codeUtils::vectorAppend(residues[0][0]->getAtoms(), residues[1][0]->getAtoms()));

    auto centerPoint      = *atoms[1]->getCoordinate();
    double maxDistance    = maxDistanceFrom(centerPoint, movingCoordinates);
    double cutoffDistance = maxDistance + constants::maxCutOff;

    auto rotationData = [&](const std::vector<Residue*>& set)
    {
        return toRotationData(Sphere {cutoffDistance, centerPoint}, dihedralResiduesMovingAtoms, set);
    };

    auto inputSets = dihedral.isBranchingLinkage ? branchedResidueSets(movingCoordinates, residues) : residues;
    return {rotationData(inputSets[0]), rotationData(inputSets[1])};
}

cds::AngleOverlap cds::wiggleWithinRangesDistanceCheck(RotatableDihedral& dihedral,
                                                       std::vector<cds::Atom*>& overlapAtomSet1,
                                                       std::vector<cds::Atom*>& overlapAtomSet2,
                                                       const DihedralAngleData& metadata, std::vector<double> angles)
{
    AngleWithMetadata initial = {cds::angle(dihedralCoordinates(dihedral)), dihedral.currentMetadata};
    std::vector<cds::AngleOverlap> results;
    for (double angle : angles)
    {
        setDihedralAngle(dihedral, {angle, metadata});
        cds::Overlap overlaps = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);

        results.push_back({
            overlaps, AngleWithMetadata {angle, metadata}
        });
    }
    // reset to initial state until this function can be made stateless
    setDihedralAngle(dihedral, initial);
    return bestOverlapResult(results);
}

// User requested gg, this prevents flipping into gt like the above would do. i.e. cb won't want a flip, gp would.
cds::AngleOverlap cds::wiggleWithinCurrentRotamer(cds::RotatableDihedral& dihedral,
                                                  std::vector<cds::Atom*>& overlapAtomSet1,
                                                  std::vector<cds::Atom*>& overlapAtomSet2, int angleIncrement)
{
    auto metadata = dihedral.currentMetadata;
    Bounds bounds = angleBounds(metadata);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    return wiggleWithinRangesDistanceCheck(dihedral, overlapAtomSet1, overlapAtomSet2, metadata,
                                           wiggleAngles(bounds, angleIncrement));
}

cds::AngleOverlap cds::wiggleUsingRotamers(const cds::DihedralCoordinates coordinates,
                                           const DihedralAngleDataVector& rotamers, int angleIncrement,
                                           const std::array<cds::DihedralRotationData, 2>& input)
{
    std::vector<AngleOverlap> results;
    for (size_t n = 0; n < rotamers.size(); n++)
    {
        const auto& rotamer = rotamers[n];
        Bounds bounds       = angleBounds(rotamer);
        AngleOverlap best   = WiggleAnglesOverlaps(coordinates, rotamer, wiggleAngles(bounds, angleIncrement), input);
        results.push_back(best);
    }

    return bestOverlapResult(results);
}
