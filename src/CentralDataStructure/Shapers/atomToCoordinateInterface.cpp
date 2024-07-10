#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <sstream>

namespace
{
    cds::Coordinate guessCoordinateOfMissingNeighbor(const cds::Atom* centralAtom, double distance)
    {
        if (centralAtom->getNeighbors().size() < 1)
        {
            std::stringstream ss;
            ss << "Error in CreateMissingCoordinateForTetrahedralAtom. centralAtom neighbors is "
               << centralAtom->getNeighbors().size() << " for " << centralAtom->getId();
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            throw std::runtime_error(ss.str());
        }
        return cds::coordinateOppositeToNeighborAverage(*centralAtom->getCoordinate(),
                                                        getCoordinatesFromAtoms(centralAtom->getNeighbors()), distance);
    }
} // namespace

// parentAtom (e.g. O of OME), childAtom (e.g. C1 of Gal1-, S1 of SO3)
void cds::moveConnectedAtomsAccordingToBondLength(cds::Atom* parentAtom, cds::Atom* childAtom)
{
    double distance      = GlycamMetadata::GetBondLengthForAtomTypes(parentAtom->getType(), childAtom->getType());
    //  Create an atom c that is will superimpose onto the a atom, bringing b atom with it.
    Coordinate c         = guessCoordinateOfMissingNeighbor(childAtom, distance);
    Coordinate cToParent = *parentAtom->getCoordinate() - c;
    // Figure out which atoms will move
    std::vector<cds::Atom*> atomsToMove;
    atomsToMove.push_back(parentAtom); // add Parent atom so search doesn't go through it.
    cdsSelections::FindConnectedAtoms(atomsToMove, childAtom);
    atomsToMove.erase(atomsToMove.begin()); // delete the parentAtom
    for (auto& atom : atomsToMove)
    {
        Coordinate* coord = atom->getCoordinate();
        *coord            = *coord + cToParent;
    }
    return;
}
