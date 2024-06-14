#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp" //FindConnectedAtoms
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Overlaps/cdsOverlaps.hpp"
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
            return std::abs(a.angle - a.metadata->default_angle_value_);
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

    cds::dihedralRotationData toRotationData(const std::vector<Atom*> movingAtoms,
                                             const std::vector<cds::Residue*> residues)
    {
        size_t atomCount = 0;
        for (auto res : residues)
        {
            atomCount += res->atomCount();
        }
        std::vector<Coordinate> coordinates;
        coordinates.reserve(atomCount);
        std::vector<std::pair<size_t, size_t>> residueAtoms;
        residueAtoms.reserve(residues.size());
        size_t currentAtom = 0;
        for (auto& res : residues)
        {
            size_t startAtom = currentAtom;
            res->insertCoordinatesInto(coordinates);
            currentAtom += res->atomCount();
            residueAtoms.push_back({startAtom, currentAtom});
        }
        std::vector<Coordinate> geometricCenters;
        geometricCenters.reserve(residues.size());
        for (size_t n = 0; n < residueAtoms.size(); n++)
        {
            auto range        = residueAtoms[n];
            Coordinate center = std::accumulate(coordinates.begin() + range.first, coordinates.begin() + range.second,
                                                Coordinate(0.0, 0.0, 0.0));
            geometricCenters.push_back(scaleBy(1.0 / (range.second - range.first), center));
        }
        const auto& firstAtoms   = residues[0]->getAtoms();
        std::vector<bool> moving = movingAtomsWithinResidue(movingAtoms, firstAtoms);
        return {coordinates, geometricCenters, residueAtoms, moving};
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
} // namespace

using cds::RotatableDihedral;

std::array<cds::dihedralRotationData, 2>
cds::dihedralRotationInputData(bool branching, const std::array<Atom*, 4>& dihedral,
                               const std::vector<Coordinate*>& movingCoordinates,
                               const std::array<std::vector<Residue*>, 2>& residues)
{
    auto dihedralResiduesMovingAtoms = movingAtomsWithinSet(
        dihedral[2], dihedral[1], codeUtils::vectorAppend(residues[0][0]->getAtoms(), residues[1][0]->getAtoms()));

    auto rotationData = [&](const std::vector<Residue*>& set)
    {
        return toRotationData(dihedralResiduesMovingAtoms, set);
    };

    auto inputSets = branching ? branchedResidueSets(movingCoordinates, residues) : residues;
    return {rotationData(inputSets[0]), rotationData(inputSets[1])};
}

double RotatableDihedral::CalculateDihedralAngle() const
{
    return cds::CalculateDihedralAngle(dihedralCoordinates());
}

int RotatableDihedral::GetNumberOfRotamers(bool likelyShapesOnly) const
{
    int count = 0;
    for (auto& entry : this->GetMetadata())
    {
        if ((entry.weight_ < 0.01) && (likelyShapesOnly)) // I'm hardcoding it, I don't care.
        {}                                                // Do nought.
        else
        {
            count++;
        }
    }
    return count;
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
    this->SetAtomsThatMove(atoms_that_move);
}

void RotatableDihedral::SetAtomsThatMove(std::vector<cds::Atom*>& atoms)
{
    coordinatesThatMove_.clear();
    for (auto& atom : atoms)
    {
        coordinatesThatMove_.push_back(atom->getCoordinate());
    }
}

std::array<Coordinate*, 4> RotatableDihedral::dihedralCoordinates() const
{
    return {atoms_[3]->getCoordinate(), atoms_[2]->getCoordinate(), atoms_[1]->getCoordinate(),
            atoms_[0]->getCoordinate()};
}

void RotatableDihedral::SetDihedralAngle(double dihedral_angle, const DihedralAngleData* metadata)
{
    if (!this->wasEverRotated_)
    {
        this->SetWasEverRotated(true);
        this->DetermineAtomsThatMove();
    }
    this->RecordPreviousState(this->CalculateDihedralAngle(), this->GetCurrentMetaData());
    auto matrix = cds::dihedralToMatrix(dihedralCoordinates(), dihedral_angle);
    matrix.rotateCoordinates(this->GetCoordinatesThatMove());
    this->SetCurrentMetaData(metadata);
}

void RotatableDihedral::SetDihedralAngleToPrevious()
{
    auto previous = this->GetPreviousState();
    this->SetDihedralAngle(previous.angle, previous.metadata);
}

double RotatableDihedral::RandomizeDihedralAngleWithinMetadataRange(const DihedralAngleData* entry)
{
    Bounds bounds = angleBounds(*entry);
    std::uniform_real_distribution<> angle_distribution(bounds.lower, bounds.upper); // define the range
    double random_angle = angle_distribution(rng);
    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/

    this->SetDihedralAngle(random_angle, entry); // THIS IS IMPORTANT!!! THIS SHOULD BE SEPARATED?!?! The two other
                                                 // functions call this one. Seems fragile.

    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/
    return random_angle;
}

void RotatableDihedral::SetRandomAngleEntryUsingMetadata()
{
    // first randomly pick one of the meta data entries. If there is only one entry, it will randomly always be it.
    std::uniform_int_distribution<> distr(0, (assigned_metadata_.size() - 1)); // define the range
    DihedralAngleData* entry = &assigned_metadata_.at(distr(rng));
    this->RandomizeDihedralAngleWithinMetadataRange(entry);
}

void RotatableDihedral::SetSpecificAngleEntryUsingMetadata(bool useRanges, long unsigned int angleEntryNumber)
{
    if (assigned_metadata_.size() <= angleEntryNumber)
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; angleEntryNumber of " << angleEntryNumber
           << " is too large as metadata.size() is " << assigned_metadata_.size() << ".\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {
        DihedralAngleData* entry = &assigned_metadata_.at(angleEntryNumber);
        if (useRanges)
        {
            this->RandomizeDihedralAngleWithinMetadataRange(entry);
        }
        else
        {
            this->SetDihedralAngle(entry->default_angle_value_, entry);
        }
    }
}

bool RotatableDihedral::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    if (dihedralName == this->GetMetadata().at(0).dihedral_angle_name_)
    {
        for (auto& metadata : this->GetMetadata())
        {
            if (metadata.rotamer_name_ == selectedRotamer)
            {
                this->SetDihedralAngle(metadata.default_angle_value_, &metadata);
                return true;
            }
        }
    }

    return false;
}

void RotatableDihedral::WiggleUsingAllRotamers(std::vector<cds::Residue*>& overlapSet1,
                                               std::vector<cds::Residue*>& overlapSet2, int angleIncrement)
{
    auto input =
        dihedralRotationInputData(isBranchingLinkage_, atoms_, coordinatesThatMove_, {overlapSet1, overlapSet2});
    auto& fixedInput                      = input[0];
    auto& movingInput                     = input[1];
    const DihedralAngleData* bestMetadata = this->GetCurrentMetaData();
    double bestDihedral                   = this->CalculateDihedralAngle();
    unsigned int lowestOverlap =
        cds::CountOverlappingAtoms({fixedInput.coordinates, fixedInput.geometricCenters, fixedInput.residueAtoms},
                                   {movingInput.coordinates, movingInput.geometricCenters, movingInput.residueAtoms});

    auto metadataVec = this->GetMetadata();
    std::vector<std::vector<AngleOverlap>> results;
    results.resize(metadataVec.size());
    for (size_t n = 0; n < metadataVec.size(); n++)
    {
        auto& metadata = metadataVec[n];
        Bounds bounds  = angleBounds(metadata);
        results[n] =
            WiggleWithinRangesDistanceCheck(fixedInput, movingInput, &metadata, wiggleAngles(bounds, angleIncrement));
    }

    for (auto& result : results)
    {
        AngleOverlap best = bestOverlapResult(result);
        if (lowestOverlap > best.overlaps)
        {
            lowestOverlap = best.overlaps;
            bestDihedral  = best.angle;
            bestMetadata  = best.metadata;
        }
    }
    this->SetDihedralAngle(bestDihedral, bestMetadata);
}

// User requested gg, this prevents flipping into gt like the above would do. i.e. cb won't want a flip, gp would.
void RotatableDihedral::WiggleWithinCurrentRotamer(std::vector<cds::Atom*>& overlapAtomSet1,
                                                   std::vector<cds::Atom*>& overlapAtomSet2, int angleIncrement)
{
    auto metadata     = GetCurrentMetaData();
    Bounds bounds     = angleBounds(*metadata);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    auto results      = this->WiggleWithinRangesDistanceCheck(overlapAtomSet1, overlapAtomSet2, metadata,
                                                              wiggleAngles(bounds, angleIncrement));
    AngleOverlap best = bestOverlapResult(results);
    this->SetDihedralAngle(best.angle, best.metadata);
}

// User requested gg, this prevents flipping into gt like the above would do.
void RotatableDihedral::WiggleWithinCurrentRotamer(std::vector<cds::Residue*>& overlapSet1,
                                                   std::vector<cds::Residue*>& overlapSet2, int angleIncrement)
{
    auto metadata = GetCurrentMetaData();
    Bounds bounds = angleBounds(*metadata);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    auto input =
        dihedralRotationInputData(isBranchingLinkage_, atoms_, coordinatesThatMove_, {overlapSet1, overlapSet2});
    auto& fixedInput  = input[0];
    auto& movingInput = input[1];
    auto results =
        this->WiggleWithinRangesDistanceCheck(fixedInput, movingInput, metadata, wiggleAngles(bounds, angleIncrement));
    AngleOverlap best = bestOverlapResult(results);
    this->SetDihedralAngle(best.angle, best.metadata);
}

std::vector<cds::AngleOverlap>
RotatableDihedral::WiggleWithinRangesDistanceCheck(std::vector<cds::Atom*>& overlapAtomSet1,
                                                   std::vector<cds::Atom*>& overlapAtomSet2,
                                                   const DihedralAngleData* metadata, std::vector<double> angles)
{
    std::vector<cds::AngleOverlap> results;
    for (double angle : angles)
    {
        this->SetDihedralAngle(angle, metadata);
        unsigned int overlaps = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);

        results.push_back({angle, metadata, overlaps});
    }
    return results;
}

std::vector<cds::AngleOverlap>
RotatableDihedral::WiggleWithinRangesDistanceCheck(const cds::dihedralRotationData& fixedInput,
                                                   const cds::dihedralRotationData& movingInput,
                                                   const DihedralAngleData* metadata, std::vector<double> angles)
{
    auto moveFirstResidueCoords = [](const RotationMatrix matrix, const dihedralRotationData& input,
                                     std::vector<Coordinate>& coordinates, std::vector<Coordinate>& centers)
    {
        auto range = input.residueAtoms[0];
        for (size_t n = range.first; n < range.second; n++)
        {
            coordinates[n] =
                input.firstResidueCoordinateMoving[n] ? matrix * input.coordinates[n] : input.coordinates[n];
        }
        Coordinate center = std::accumulate(coordinates.begin() + range.first, coordinates.begin() + range.second,
                                            Coordinate(0.0, 0.0, 0.0));
        centers[0]        = scaleBy(1.0 / (range.second - range.first), center);
    };
    const auto dihedral                       = dihedralCoordinates();
    // copies of input vectors used for updates during looping
    std::vector<Coordinate> fixedCoordinates  = fixedInput.coordinates;
    std::vector<Coordinate> fixedCenters      = fixedInput.geometricCenters;
    std::vector<Coordinate> movingCoordinates = movingInput.coordinates;
    std::vector<Coordinate> movingCenters     = movingInput.geometricCenters;
    std::vector<cds::AngleOverlap> results;
    for (double angle : angles)
    {
        auto matrix = cds::dihedralToMatrix(dihedral, angle);
        moveFirstResidueCoords(matrix, fixedInput, fixedCoordinates, fixedCenters);
        moveFirstResidueCoords(matrix, movingInput, movingCoordinates, movingCenters);
        for (size_t n = movingInput.residueAtoms[0].second; n < movingCoordinates.size(); n++)
        {
            movingCoordinates[n] = matrix * movingInput.coordinates[n];
        }
        for (size_t n = 1; n < movingCenters.size(); n++)
        {
            movingCenters[n] = matrix * movingInput.geometricCenters[n];
        }
        unsigned int overlaps = CountOverlappingAtoms({fixedCoordinates, fixedCenters, fixedInput.residueAtoms},
                                                      {movingCoordinates, movingCenters, movingInput.residueAtoms});

        results.push_back({angle, metadata, overlaps});
    }
    return results;
}

std::string RotatableDihedral::Print() const
{
    std::stringstream ss;
    ss << atoms_[0]->getName() << ", " << atoms_[1]->getName() << ", " << atoms_[2]->getName() << ", "
       << atoms_[3]->getName() << ": " << this->CalculateDihedralAngle() << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
}
