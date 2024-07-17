#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp" //FindConnectedAtoms
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>

namespace
{
    struct Bounds
    {
        double lower;
        double upper;
    };

    struct Dihedral
    {
        double angle;
        unsigned int overlaps;
    };

    Bounds angleBounds(const DihedralAngleData& metadata)
    {
        double defaultAngle = metadata.default_angle_value_;
        return {defaultAngle - metadata.lower_deviation_, defaultAngle + metadata.upper_deviation_};
    }

    cds::AngleWithMetadata randomDihedralAngleWithinMetadataRange(const DihedralAngleData* entry)
    {
        Bounds bounds = angleBounds(*entry);
        std::uniform_real_distribution<> angle_distribution(bounds.lower, bounds.upper); // define the range
        double random_angle = angle_distribution(rng);
        return {random_angle, entry};
    }

    std::vector<double> wiggleAngles(Bounds bounds, int approximateIncrement)
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
            return std::abs(a.angle.value - a.angle.metadata->default_angle_value_);
        };
        int best = 0;
        for (size_t n = 1; n < results.size(); n++)
        {
            auto a = results[n];
            auto b = results[best];
            if ((a.overlaps < b.overlaps) ||
                ((a.overlaps == b.overlaps) && differenceFromDefault(a) < differenceFromDefault(b)))
            {
                best = n;
            }
        }
        return results[best];
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

    std::vector<bool> movingAtomsWithinResidue(const std::vector<Atom*> movingAtoms,
                                               const std::vector<Atom*> residueAtoms)
    {
        std::vector<bool> result;
        for (auto& a : residueAtoms)
        {
            result.push_back(codeUtils::contains(movingAtoms, a));
        }
        return result;
    }

    cds::dihedralRotationData toRotationData(const cds::Sphere bounds, const std::vector<cds::Atom*> movingAtoms,
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

    std::array<std::vector<cds::Residue*>, 2> branchedResidueSets(const std::vector<Coordinate*>& movingCoordinates,
                                                                  const std::array<std::vector<Residue*>, 2>& residues)
    {
        auto residueContainsMovingAtom = [&](Residue* res)
        {
            return codeUtils::contains(movingCoordinates, res->getAtoms()[0]->getCoordinate());
        };
        auto& setA                    = residues[0];
        auto& setB                    = residues[1];
        std::vector<Residue*> newSetA = setB;
        std::vector<Residue*> newSetB {setA[0]};
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

    void moveFirstResidueCoords(const cds::RotationMatrix matrix, const cds::dihedralRotationData& input,
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

    cds::AngleOverlap WiggleAnglesOverlaps(const cds::DihedralCoordinates dihedral,
                                           const cds::dihedralRotationData& fixedInput,
                                           const cds::dihedralRotationData& movingInput,
                                           const DihedralAngleData* metadata, std::vector<double> angles)
    {
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
            unsigned int overlaps =
                cds::CountOverlappingAtoms({fixedCoordinates, fixedSpheres, fixedInput.residueAtoms},
                                           {movingCoordinates, movingSpheres, movingInput.residueAtoms});

            results.push_back({
                overlaps, cds::AngleWithMetadata {angle, metadata}
            });
        }
        return bestOverlapResult(results);
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
} // namespace

using cds::RotatableDihedral;

std::array<cds::dihedralRotationData, 2>
cds::dihedralRotationInputData(bool branching, const std::array<Atom*, 4>& dihedral,
                               const std::vector<Coordinate*>& movingCoordinates,
                               const std::array<std::vector<Residue*>, 2>& residues)
{
    auto dihedralResiduesMovingAtoms = movingAtomsWithinSet(
        dihedral[2], dihedral[1], codeUtils::vectorAppend(residues[0][0]->getAtoms(), residues[1][0]->getAtoms()));

    auto centerPoint      = *dihedral[1]->getCoordinate();
    double maxDistance    = maxDistanceFrom(centerPoint, movingCoordinates);
    double cutoffDistance = maxDistance + constants::maxCutOff;

    auto rotationData = [&](const std::vector<Residue*>& set)
    {
        return toRotationData(Sphere {cutoffDistance, centerPoint}, dihedralResiduesMovingAtoms, set);
    };

    auto inputSets = branching ? branchedResidueSets(movingCoordinates, residues) : residues;
    return {rotationData(inputSets[0]), rotationData(inputSets[1])};
}

int RotatableDihedral::GetNumberOfRotamers(bool likelyShapesOnly) const
{
    return likelyShapesOnly ? GetLikelyMetadata().size() : GetMetadata().size();
}

DihedralAngleDataVector RotatableDihedral::GetLikelyMetadata() const
{
    DihedralAngleDataVector returningMetadata;
    returningMetadata.reserve(assigned_metadata_.size());
    for (auto& entry : assigned_metadata_)
    {
        if (entry.weight_ >= 0.01) // HARDCODE EVERYTHING.
        {
            returningMetadata.push_back(entry);
        }
    }
    return returningMetadata;
}

std::string RotatableDihedral::GetName() const
{
    if (this->GetLikelyMetadata().empty())
    {
        return "Boo";
    }
    else
    {
        return this->GetLikelyMetadata().at(0).dihedral_angle_name_;
    }
}

void RotatableDihedral::DetermineAtomsThatMove()
{
    // In keeping with giving residues as GlcNAc1-4Gal, and wanting the moving atoms to be in the opposite direction
    // (defaults):
    std::vector<cds::Atom*> atoms_that_move;
    atoms_that_move.push_back(atoms_[2]);
    cdsSelections::FindConnectedAtoms(atoms_that_move, atoms_[1]);
    atoms_that_move.erase(atoms_that_move.begin());
    coordinatesThatMove_ = getCoordinatesFromAtoms(atoms_that_move);
}

cds::DihedralCoordinates RotatableDihedral::dihedralCoordinates() const
{
    std::array<Coordinate*, 4> coords = {atoms_[3]->getCoordinate(), atoms_[2]->getCoordinate(),
                                         atoms_[1]->getCoordinate(), atoms_[0]->getCoordinate()};
    return DihedralCoordinates {*coords[0], *coords[1], *coords[2], *coords[3]};
}

void RotatableDihedral::SetDihedralAngle(AngleWithMetadata target)
{
    this->RecordPreviousState({cds::angle(dihedralCoordinates()), GetCurrentMetaData()});
    auto matrix = rotationTo(dihedralCoordinates(), constants::degree2Radian(target.value));
    matrix.rotateCoordinates(this->GetCoordinatesThatMove());
    this->SetCurrentMetaData(target.metadata);
}

cds::AngleWithMetadata RotatableDihedral::RandomAngleEntryUsingMetadata()
{
    // first randomly pick one of the meta data entries. If there is only one entry, it will randomly always be it.
    std::uniform_int_distribution<> distr(0, (assigned_metadata_.size() - 1)); // define the range
    DihedralAngleData* entry = &assigned_metadata_.at(distr(rng));
    return randomDihedralAngleWithinMetadataRange(entry);
}

cds::AngleWithMetadata RotatableDihedral::SpecificAngleEntryUsingMetadata(bool useRanges,
                                                                          long unsigned int angleEntryNumber)
{
    if (assigned_metadata_.size() <= angleEntryNumber)
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; angleEntryNumber of " << angleEntryNumber
           << " is too large as metadata.size() is " << assigned_metadata_.size() << ".\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    DihedralAngleData* entry = &assigned_metadata_.at(angleEntryNumber);
    return useRanges ? randomDihedralAngleWithinMetadataRange(entry)
                     : AngleWithMetadata {entry->default_angle_value_, entry};
}

bool RotatableDihedral::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    if (dihedralName == this->GetMetadata().at(0).dihedral_angle_name_)
    {
        for (auto& metadata : this->GetMetadata())
        {
            if (metadata.rotamer_name_ == selectedRotamer)
            {
                this->SetDihedralAngle({metadata.default_angle_value_, &metadata});
                return true;
            }
        }
    }

    return false;
}

cds::AngleOverlap RotatableDihedral::WiggleWithinRangesDistanceCheck(std::vector<cds::Atom*>& overlapAtomSet1,
                                                                     std::vector<cds::Atom*>& overlapAtomSet2,
                                                                     const DihedralAngleData* metadata,
                                                                     std::vector<double> angles)
{
    AngleWithMetadata initial = {cds::angle(dihedralCoordinates()), GetCurrentMetaData()};
    std::vector<cds::AngleOverlap> results;
    for (double angle : angles)
    {
        SetDihedralAngle({angle, metadata});
        unsigned int overlaps = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);

        results.push_back({
            overlaps, AngleWithMetadata {angle, metadata}
        });
    }
    // reset to initial state until this function can be made stateless
    SetDihedralAngle(initial);
    return bestOverlapResult(results);
}

// User requested gg, this prevents flipping into gt like the above would do. i.e. cb won't want a flip, gp would.
cds::AngleOverlap RotatableDihedral::WiggleWithinCurrentRotamer(std::vector<cds::Atom*>& overlapAtomSet1,
                                                                std::vector<cds::Atom*>& overlapAtomSet2,
                                                                int angleIncrement)
{
    auto metadata = GetCurrentMetaData();
    Bounds bounds = angleBounds(*metadata);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    return WiggleWithinRangesDistanceCheck(overlapAtomSet1, overlapAtomSet2, metadata,
                                           wiggleAngles(bounds, angleIncrement));
}

cds::AngleOverlap RotatableDihedral::WiggleUsingRotamers(const DihedralAngleDataVector& rotamers, int angleIncrement,
                                                         const std::array<std::vector<cds::Residue*>, 2>& residues)
{
    auto input        = dihedralRotationInputData(isBranchingLinkage_, atoms_, coordinatesThatMove_, residues);
    auto& fixedInput  = input[0];
    auto& movingInput = input[1];
    std::vector<AngleOverlap> results;
    for (size_t n = 0; n < rotamers.size(); n++)
    {
        const auto& rotamer = rotamers[n];
        Bounds bounds       = angleBounds(rotamer);
        AngleOverlap best   = WiggleAnglesOverlaps(dihedralCoordinates(), fixedInput, movingInput, &rotamer,
                                                   wiggleAngles(bounds, angleIncrement));
        results.push_back(best);
    }

    return bestOverlapResult(results);
}

std::string RotatableDihedral::Print() const
{
    std::stringstream ss;
    ss << atoms_[0]->getName() << ", " << atoms_[1]->getName() << ", " << atoms_[2]->getName() << ", "
       << atoms_[3]->getName() << ": " << cds::angle(dihedralCoordinates()) << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
}
