#include "includes/CentralDataStructure/Shapers/psiAngleHydrogen.hpp"

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomCoordinateInterface.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <vector>
#include <memory>

cds::Atom* cds::findHydrogenForPsiAngle(const Atom* atom)
{
    for (auto& neighbor : atom->getNeighbors())
    {
        if (neighbor->getName().at(0) == 'H')
        {
            return neighbor;
        }
    }
    return nullptr;
}

void cds::createHydrogenForPsiAngles(Residue* residue, std::vector<DihedralAtoms>& dihedralAtoms,
                                     const DihedralAngleMetadata& metadata)
{
    for (size_t n = 0; n < dihedralAtoms.size(); n++)
    {
        for (auto& entry : metadata[n])
        {
            if (entry.dihedral_angle_name_ == "Psi" && entry.atom4_.at(0) == 'H')
            { // If it's a psi angle and is supposed to be defined by a H...
                Atom* atom     = dihedralAtoms[n].atoms[2];
                Atom* hydrogen = findHydrogenForPsiAngle(atom);
                if (hydrogen != nullptr)
                {
                    dihedralAtoms[n].atoms[3] = hydrogen;
                }
                else
                {
                    Coordinate newCoord = coordinateOppositeToNeighborAverage(
                        *atom->getCoordinate(), getCoordinatesFromAtoms(atom->getNeighbors()), 1.0);
                    std::unique_ptr<Atom> newAtom = std::make_unique<Atom>("HHH", newCoord);
                    addBond(atom, newAtom.get());
                    residue->addAtom(std::move(newAtom));
                }
            }
        }
    }
}
