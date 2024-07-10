#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <sstream>

cds::Coordinate cds::GuessMissingCoordinateForAtom(cds::Atom* centralAtom, const double distance)
{
    if (centralAtom->getNeighbors().size() < 1)
    {
        std::stringstream ss;
        ss << "Error in CreateMissingCoordinateForTetrahedralAtom. centralAtom neighbors is "
           << centralAtom->getNeighbors().size() << " for " << centralAtom->getId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    std::vector<Coordinate*> threeNeighborCoords;
    std::vector<Atom*> neighbors = centralAtom->getNeighbors();
    for (auto& neighbor : neighbors)
    {
        //        std::cout << neighbor->getName() << ", " << neighbor->getCoordinate()->ToString() << "\n";
        threeNeighborCoords.push_back(neighbor->getCoordinate());
    }
    return cds::CreateCoordinateForCenterAwayFromNeighbors(*centralAtom->getCoordinate(), threeNeighborCoords,
                                                           distance);
}

// parentAtom (e.g. O of OME), childAtom (e.g. C1 of Gal1-, S1 of SO3)
void cds::FindAtomsToMoveSetDistance(cds::Atom* parentAtom, cds::Atom* childAtom)
{ // Figure out distance
    //    std::stringstream ss;
    //    ss << "parent is " << parentAtom->getName() << "_" << parentAtom->getIndex() << " "
    //       << parentAtom->getCoordinate()->ToString() << "\n";
    //    ss << "child is " << childAtom->getName() << "_" << childAtom->getIndex() << " "
    //       << childAtom->getCoordinate()->ToString() << "\n";
    double distance = GlycamMetadata::GetBondLengthForAtomTypes(parentAtom->getType(), childAtom->getType());
    //  Create an atom c that is will superimpose onto the a atom, bringing b atom with it.
    Coordinate c    = cds::GuessMissingCoordinateForAtom(childAtom, distance);
    // std::cout << "New tetraAtom for child is: " << c.ToString() << "\n";
    Coordinate cToParent(parentAtom->getCoordinate()->GetX() - c.GetX(), parentAtom->getCoordinate()->GetY() - c.GetY(),
                         parentAtom->getCoordinate()->GetZ() - c.GetZ());
    // Figure out which atoms will move
    std::vector<cds::Atom*> atomsToRotate;
    atomsToRotate.push_back(parentAtom); // add Parent atom so search doesn't go through it.
    cdsSelections::FindConnectedAtoms(atomsToRotate, childAtom);
    atomsToRotate.erase(atomsToRotate.begin()); // delete the parentAtom
    for (auto& atom : atomsToRotate)
    {
        Coordinate* coord = atom->getCoordinate();
        *coord            = *coord + cToParent;
        //        ss << "Moved " << atom->getName() << "_" << atom->getIndex() << " to new position:\n";
        //        atom->getCoordinate()->Print(ss);
        //        ss << "\n";
    }
    //    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return;
}
