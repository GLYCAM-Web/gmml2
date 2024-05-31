#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp" // getCoordinatesFromAtoms
#include "includes/CentralDataStructure/Shapers/shapers.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp" //FindConnectedAtoms
#include "includes/CentralDataStructure/Measurements/measurements.hpp"

#include <sstream>

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

namespace
{
    Bounds angleBounds(const DihedralAngleData& metadata)
    {
        double defaultAngle = metadata.default_angle_value_;
        return {defaultAngle - metadata.lower_deviation_, defaultAngle + metadata.upper_deviation_};
    }
} // namespace

using cds::RotatableDihedral;

double RotatableDihedral::CalculateDihedralAngle(const std::string type) const
{
    if (type == "glycamReport")
    {
        if ((this->GetName() == "Phi") && (atoms_[0]->getName() == "C1"))
        { // ToDO there will be an issue with C1-O5 linkages unless atom knows which residue it is in.
            // Coordinate* o5Coord = atoms_[0]->GetResidue()->GetAtom("O5")->getCoordinate();
            Coordinate* o5Coord = cdsSelections::getNeighborNamed(atoms_[0], "O5")->getCoordinate();
            return cds::CalculateDihedralAngle(
                {o5Coord, atoms_[1]->getCoordinate(), atoms_[2]->getCoordinate(), atoms_[3]->getCoordinate()});
        }
    }
    return cds::CalculateDihedralAngle({atoms_[0]->getCoordinate(), atoms_[1]->getCoordinate(),
                                        atoms_[2]->getCoordinate(), atoms_[3]->getCoordinate()});
}

int RotatableDihedral::GetNumberOfRotamers(bool likelyShapesOnly) const
{
    if (this->GetMetadata().empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Error in RotatableDihedral::GetNumberOfRotamers; no metadata has been set.\n");
        return 0;
    }
    else
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
}

DihedralAngleDataVector RotatableDihedral::GetLikelyMetadata() const
{
    DihedralAngleDataVector returningMetadata;
    for (auto& entry : this->GetMetadata())
    {
        if (entry.weight_ >= 0.01) // HARDCODE EVERYTHING.
        {
            returningMetadata.push_back(entry);
        }
    }
    return returningMetadata;
}

std::vector<double> RotatableDihedral::GetAllPossibleAngleValues(const int interval) const
{
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::GetAllPossibleAngleValues; no metadata has been set for:\n"
           << atoms_[0]->getId() << " " << atoms_[1]->getId() << " " << atoms_[2]->getId() << " " << atoms_[3]->getId()
           << "\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {
        std::vector<double> allPossibleAngleValues;
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_entries = this->GetMetadata();
        for (auto& metadata : metadata_entries)
        {
            Bounds bounds = angleBounds(metadata);
            for (double angle = bounds.lower; angle <= bounds.upper; angle += interval)
            {
                allPossibleAngleValues.push_back(angle);
            }
        }
        return allPossibleAngleValues;
    }
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
    int moving = GetIsAtomsThatMoveReversed() ? 2 : 1;
    int fixed  = GetIsAtomsThatMoveReversed() ? 1 : 2;
    atoms_that_move.push_back(atoms_[moving]);
    cdsSelections::FindConnectedAtoms(atoms_that_move, atoms_[fixed]);
    this->SetAtomsThatMove(atoms_that_move);
}

void RotatableDihedral::SetAtomsThatMove(std::vector<cds::Atom*>& atoms)
{
    coordinatesThatMove_.clear();
    this->AddExtraAtomsThatMove(atoms);
}

void RotatableDihedral::AddExtraAtomsThatMove(std::vector<cds::Atom*>& extraAtoms)
{
    for (auto& atom : extraAtoms)
    {
        coordinatesThatMove_.push_back(atom->getCoordinate());
    }
}

void RotatableDihedral::SetDihedralAngle(const double dihedral_angle)
{
    if (!this->wasEverRotated_)
    {
        this->SetWasEverRotated(true);
        this->DetermineAtomsThatMove();
    }
    Coordinate* a1 = atoms_[0]->getCoordinate();
    Coordinate* a2 = atoms_[1]->getCoordinate();
    Coordinate* a3 = atoms_[2]->getCoordinate();
    Coordinate* a4 = atoms_[3]->getCoordinate();
    if (this->GetIsAtomsThatMoveReversed())
    { // A function that was agnostic of this would be great.
        a1 = atoms_[3]->getCoordinate();
        a2 = atoms_[2]->getCoordinate();
        a3 = atoms_[1]->getCoordinate();
        a4 = atoms_[0]->getCoordinate();
    }
    this->RecordPreviousDihedralAngle(this->CalculateDihedralAngle());
    cds::SetDihedralAngle(a1, a2, a3, a4, dihedral_angle, this->GetCoordinatesThatMove());
}

void RotatableDihedral::SetDihedralAngleToPrevious()
{
    this->SetDihedralAngle(this->GetPreviousDihedralAngle());
}

double RotatableDihedral::RandomizeDihedralAngle()
{
    return RotatableDihedral::RandomizeDihedralAngleWithinRange(0.0, 360.0);
}

double RotatableDihedral::RandomizeDihedralAngleWithinRange(double min, double max)
{
    std::uniform_real_distribution<> angle_distribution(min, max); // define the range
    double random_angle = angle_distribution(rng);
    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/

    this->SetDihedralAngle(random_angle); // THIS IS IMPORTANT!!! THIS SHOULD BE SEPARATED?!?! The two other functions
                                          // call this one. Seems fragile.

    /*******************************************/
    /*               IMPORTANT                 */
    /*******************************************/
    return random_angle;
}

void RotatableDihedral::AddMetadata(DihedralAngleData metadata)
{
    assigned_metadata_.push_back(metadata);
}

void RotatableDihedral::SetRandomAngleEntryUsingMetadata(bool useRanges)
{
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetRandomAngleEntryUsingMetadata; no metadata has been set for:.\n"
           << atoms_[0]->getId() << atoms_[1]->getId() << atoms_[2]->getId() << atoms_[3]->getId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    { // first randomly pick one of the meta data entries. If there is only one entry, it will randomly always be it.
        std::uniform_int_distribution<> distr(0, (assigned_metadata_.size() - 1)); // define the range
        DihedralAngleData& entry = assigned_metadata_.at(distr(rng));
        this->SetCurrentMetaData(entry);
        if (useRanges)
        {
            Bounds bounds = angleBounds(entry);
            this->RandomizeDihedralAngleWithinRange(bounds.lower, bounds.upper);
        }
        else
        { // Just set it to the entries default value
            this->RandomizeDihedralAngleWithinRange(entry.default_angle_value_, entry.default_angle_value_);
        }
    }
}

void RotatableDihedral::SetSpecificAngleEntryUsingMetadata(bool useRanges, long unsigned int angleEntryNumber)
{
    if (assigned_metadata_.empty())
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; no metadata has been set for:\n"
           << atoms_[0]->getId() << " " << atoms_[1]->getId() << " " << atoms_[2]->getId() << " " << atoms_[3]->getId()
           << "\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else if (assigned_metadata_.size() <= angleEntryNumber)
    {
        std::stringstream ss;
        ss << "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; angleEntryNumber of " << angleEntryNumber
           << " is too large as metadata.size() is " << assigned_metadata_.size() << ".\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    else
    {
        DihedralAngleData& entry = assigned_metadata_.at(angleEntryNumber);
        this->SetCurrentMetaData(entry);
        if (useRanges)
        {
            Bounds bounds = angleBounds(entry);
            this->RandomizeDihedralAngleWithinRange(bounds.lower, bounds.upper);
        }
        else
        {
            this->SetDihedralAngle(entry.default_angle_value_);
        }
    }
}

bool RotatableDihedral::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    if (assigned_metadata_.empty())
    {
        std::string message =
            "Error in RotatableDihedral::SetSpecificAngleUsingMetadata; no metadata has been set for:\n" +
            atoms_[0]->getId() + " " + atoms_[1]->getId() + " " + atoms_[2]->getId() + " " + atoms_[3]->getId();
        message += ".\nA rotatable dihedral named: " + dihedralName +
                   " was requested to be set to the rotamer named: " + selectedRotamer;
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    if (dihedralName == this->GetMetadata().at(0).dihedral_angle_name_)
    {
        for (auto& metadata : this->GetMetadata())
        {
            if (metadata.rotamer_name_ == selectedRotamer)
            {
                this->SetDihedralAngle(metadata.default_angle_value_);
                this->SetCurrentMetaData(metadata);
                return true;
            }
        }
    }

    return false;
}

void RotatableDihedral::WiggleUsingAllRotamers(std::vector<cds::Residue*>& overlapSet1,
                                               std::vector<cds::Residue*>& overlapSet2, const int& angleIncrement)
{

    double bestDihedral        = this->CalculateDihedralAngle();
    auto overlapInput          = cds::toResidueAtomOverlapInput(overlapSet1, overlapSet2);
    unsigned int lowestOverlap = cds::CountOverlappingAtoms(overlapInput.first, overlapInput.second);
    for (auto& metadata : this->GetMetadata())
    {
        Bounds bounds = angleBounds(metadata);
        unsigned int newOverlap =
            this->WiggleWithinRangesDistanceCheck(overlapInput, angleIncrement, bounds.lower, bounds.upper);
        if (lowestOverlap >= (newOverlap + 1))
        {
            lowestOverlap = newOverlap;
            bestDihedral  = this->CalculateDihedralAngle();
            this->SetCurrentMetaData(metadata);
        }
        else
        {
            this->SetDihedralAngle(bestDihedral);
        }
    }
}

void RotatableDihedral::WiggleUsingAllRotamers(std::vector<cds::Atom*>& overlapAtomSet1,
                                               std::vector<cds::Atom*>& overlapAtomSet2, const int& angleIncrement)
{
    double bestDihedral        = this->CalculateDihedralAngle();
    unsigned int lowestOverlap = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);
    ;
    for (auto& metadata : this->GetMetadata())
    {
        Bounds bounds           = angleBounds(metadata);
        unsigned int newOverlap = this->WiggleWithinRangesDistanceCheck(overlapAtomSet1, overlapAtomSet2,
                                                                        angleIncrement, bounds.lower, bounds.upper);
        if (lowestOverlap >= (newOverlap + 1))
        {
            lowestOverlap = newOverlap;
            bestDihedral  = this->CalculateDihedralAngle();
            this->SetCurrentMetaData(metadata);
        }
        else
        {
            this->SetDihedralAngle(bestDihedral);
        }
    }
}

// User requested gg, this prevents flipping into gt like the above would do. i.e. cb won't want a flip, gp would.
void RotatableDihedral::WiggleWithinCurrentRotamer(std::vector<cds::Atom*>& overlapAtomSet1,
                                                   std::vector<cds::Atom*>& overlapAtomSet2, const int& angleIncrement)
{
    Bounds bounds = angleBounds(*GetCurrentMetaData());
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    this->WiggleWithinRangesDistanceCheck(overlapAtomSet1, overlapAtomSet2, angleIncrement, bounds.lower, bounds.upper);
}

// User requested gg, this prevents flipping into gt like the above would do.
void RotatableDihedral::WiggleWithinCurrentRotamer(std::vector<cds::Residue*>& overlapResidueSet1,
                                                   std::vector<cds::Residue*>& overlapResidueSet2,
                                                   const int& angleIncrement)
{
    Bounds bounds = angleBounds(*GetCurrentMetaData());
    // Set to lowest deviation, work through to highest. Set best value and return it for reference.
    auto input    = toResidueAtomOverlapInput(overlapResidueSet1, overlapResidueSet2);
    this->WiggleWithinRangesDistanceCheck(input, angleIncrement, bounds.lower, bounds.upper);
}

unsigned int RotatableDihedral::WiggleWithinRangesDistanceCheck(std::vector<cds::Atom*>& overlapAtomSet1,
                                                                std::vector<cds::Atom*>& overlapAtomSet2,
                                                                const int& angleIncrement, const double& lowerBound,
                                                                const double& upperBound)
{
    this->SetDihedralAngle(lowerBound);
    double currentDihedral      = lowerBound;
    double bestDihedral         = lowerBound;
    unsigned int currentOverlap = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);
    unsigned int lowestOverlap  = currentOverlap;
    double defaultAngle         = GetCurrentMetaData()->default_angle_value_;
    while (currentDihedral < upperBound)
    {
        currentDihedral += angleIncrement; // increment
        this->SetDihedralAngle(currentDihedral);
        currentOverlap = cds::CountOverlappingAtoms(overlapAtomSet1, overlapAtomSet2);
        if (currentOverlap < lowestOverlap)
        {
            lowestOverlap = currentOverlap;
            bestDihedral  = currentDihedral;
        }
        // Prefer angles closer to default if overlap is the same.
        else if ((lowestOverlap == currentOverlap) &&
                 (std::abs(defaultAngle - bestDihedral) > std::abs(defaultAngle - currentDihedral)))
        {
            bestDihedral = currentDihedral;
        }
    }
    this->SetDihedralAngle(bestDihedral);
    return lowestOverlap;
}

unsigned int RotatableDihedral::WiggleWithinRangesDistanceCheck(cds::ResidueAtomOverlapInputPair& overlapInput,
                                                                const int& angleIncrement, const double& lowerBound,
                                                                const double& upperBound)
{
    double currentDihedral = lowerBound;
    double bestDihedral    = lowerBound;
    this->SetDihedralAngle(lowerBound);
    cds::setGeometricCenters(overlapInput);
    double defaultAngle = GetCurrentMetaData()->default_angle_value_;

    unsigned int currentOverlap = cds::CountOverlappingAtoms(overlapInput.first, overlapInput.second);
    unsigned int lowestOverlap  = currentOverlap;
    if (this->GetCurrentMetaData() == nullptr) // you don't need these checks if you have RAII OLIVER
    {
        std::string message = "Error: current metadata not set for RotatableDihedral " + this->GetName();
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    while (currentDihedral < upperBound)
    {
        currentDihedral += angleIncrement; // increment
        this->SetDihedralAngle(currentDihedral);
        cds::setGeometricCenters(overlapInput);

        currentOverlap = cds::CountOverlappingAtoms(overlapInput.first, overlapInput.second);
        if (currentOverlap < lowestOverlap)
        {
            lowestOverlap = currentOverlap;
            bestDihedral  = currentDihedral;
        }
        // Prefer angles closer to default if overlap is the same.
        else if ((lowestOverlap == currentOverlap) &&
                 (std::abs(defaultAngle - bestDihedral) > std::abs(defaultAngle - currentDihedral)))
        {
            bestDihedral = currentDihedral;
        }
    }
    this->SetDihedralAngle(bestDihedral);
    return lowestOverlap;
}

bool RotatableDihedral::IsThereHydrogenForPsiAngle()
{
    for (auto& neighbor : atoms_[2]->getNeighbors())
    {
        if (neighbor->getName().at(0) == 'H')
        {
            atoms_[3] = neighbor;
            return true;
        }
    }
    return false;
}

// Private

std::unique_ptr<cds::Atom> RotatableDihedral::CreateHydrogenAtomForPsiAngle()
{
    std::vector<Coordinate*> neighborsCoords = cds::getCoordinatesFromAtoms(atoms_[2]->getNeighbors());
    Coordinate newCoord = cds::CreateCoordinateForCenterAwayFromNeighbors(*atoms_[2]->getCoordinate(), neighborsCoords);
    std::unique_ptr<cds::Atom> newAtom = std::make_unique<cds::Atom>("HHH", newCoord);
    atoms_[2]->addBond(newAtom.get());
    return newAtom;
}

std::string RotatableDihedral::Print() const
{
    std::stringstream ss;
    ss << atoms_[0]->getName() << ", " << atoms_[1]->getName() << ", " << atoms_[2]->getName() << ", "
       << atoms_[3]->getName() << ": " << this->CalculateDihedralAngle() << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
}
