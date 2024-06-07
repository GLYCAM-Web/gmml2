#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp" //FindConnectedAtoms
#include "includes/CentralDataStructure/Measurements/measurements.hpp"

#include <algorithm>
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

    cds::AngleOverlap bestOverlapResult(const std::vector<cds::AngleOverlap>& results, double defaultAngle)
    {
        int best = 0;
        for (size_t n = 1; n < results.size(); n++)
        {
            auto a = results[n];
            auto b = results[best];
            if ((a.overlaps < b.overlaps) ||
                ((a.overlaps == b.overlaps) && std::abs(a.angle - defaultAngle) < std::abs(b.angle - defaultAngle)))
            {
                best = n;
            }
        }
        return results[best];
    }

} // namespace

using cds::RotatableDihedral;

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
    bool reversed                         = overlapSet2[0]->contains(atoms_[3]);
    const DihedralAngleData* bestMetadata = this->GetCurrentMetaData();
    double bestDihedral                   = this->CalculateDihedralAngle();
    auto overlapInput                     = reversed ? cds::toResidueAtomOverlapInput(overlapSet2, overlapSet1, true)
                                                     : cds::toResidueAtomOverlapInput(overlapSet1, overlapSet2, true);
    unsigned int lowestOverlap            = cds::CountOverlappingAtoms(overlapInput.first, overlapInput.second);

    for (auto& metadata : this->GetMetadata())
    {
        Bounds bounds = angleBounds(metadata);
        auto results =
            this->WiggleWithinRangesDistanceCheck(overlapInput, &metadata, wiggleAngles(bounds, angleIncrement));
        AngleOverlap best = bestOverlapResult(results, metadata.default_angle_value_);
        if (lowestOverlap > best.overlaps)
        {
            lowestOverlap = best.overlaps;
            bestDihedral  = best.angle;
            bestMetadata  = &metadata;
        }
    }
    this->SetDihedralAngle(bestDihedral, bestMetadata);
}

void RotatableDihedral::WiggleUsingAllRotamers(std::vector<cds::Atom*>& overlapAtomSet1,
                                               std::vector<cds::Atom*>& overlapAtomSet2, int angleIncrement)
{
    const DihedralAngleData* bestMetadata = this->GetCurrentMetaData();
    double bestDihedral                   = this->CalculateDihedralAngle();
    unsigned int lowestOverlap            = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);

    for (auto& metadata : this->GetMetadata())
    {
        Bounds bounds     = angleBounds(metadata);
        auto results      = this->WiggleWithinRangesDistanceCheck(overlapAtomSet1, overlapAtomSet2, &metadata,
                                                                  wiggleAngles(bounds, angleIncrement));
        AngleOverlap best = bestOverlapResult(results, metadata.default_angle_value_);
        if (lowestOverlap > best.overlaps)
        {
            lowestOverlap = best.overlaps;
            bestDihedral  = best.angle;
            bestMetadata  = &metadata;
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
    AngleOverlap best = bestOverlapResult(results, metadata->default_angle_value_);
    this->SetDihedralAngle(best.angle, metadata);
}

// User requested gg, this prevents flipping into gt like the above would do.
void RotatableDihedral::WiggleWithinCurrentRotamer(std::vector<cds::Residue*>& overlapResidueSet1,
                                                   std::vector<cds::Residue*>& overlapResidueSet2, int angleIncrement)
{
    auto metadata     = GetCurrentMetaData();
    Bounds bounds     = angleBounds(*metadata);
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    auto input        = toResidueAtomOverlapInput(overlapResidueSet1, overlapResidueSet2, false);
    auto results      = this->WiggleWithinRangesDistanceCheck(input, metadata, wiggleAngles(bounds, angleIncrement));
    AngleOverlap best = bestOverlapResult(results, metadata->default_angle_value_);
    this->SetDihedralAngle(best.angle, metadata);
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

        results.push_back({angle, overlaps});
    }
    return results;
}

std::vector<cds::AngleOverlap>
RotatableDihedral::WiggleWithinRangesDistanceCheck(cds::ResidueAtomOverlapInputPair& overlapInput,
                                                   const DihedralAngleData* metadata, std::vector<double> angles)
{
    std::vector<cds::AngleOverlap> results;
    for (double angle : angles)
    {
        this->SetDihedralAngle(angle, metadata);
        cds::setGeometricCenters(overlapInput);
        unsigned int overlaps = cds::CountOverlappingAtoms(overlapInput.first, overlapInput.second);

        results.push_back({angle, overlaps});
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
