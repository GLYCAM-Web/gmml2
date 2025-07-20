#include "includes/CentralDataStructure/Shapers/psiAngleHydrogen.hpp"

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <memory>
#include <vector>

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

void cds::createHydrogenForPsiAngles(
    const GlycamMetadata::DihedralAngleDataTable& metadataTable,
    Residue* residue,
    std::vector<DihedralAtoms>& dihedralAtoms,
    const std::vector<std::vector<size_t>>& metadataIndices)
{
    for (size_t n = 0; n < dihedralAtoms.size(); n++)
    {
        for (auto& entry : metadataIndices[n])
        {
            if (metadataTable.entries[entry].dihedral_angle_name_ == "Psi" &&
                metadataTable.entries[entry].atom4_.at(0) == 'H')
            { // If it's a psi angle and is supposed to be defined by a H...
                Atom* atom = dihedralAtoms[n][2];
                Atom* hydrogen = findHydrogenForPsiAngle(atom);
                if (hydrogen != nullptr)
                {
                    dihedralAtoms[n][3] = hydrogen;
                }
                else
                {
                    Coordinate newCoord = coordinateOppositeToNeighborAverage(
                        atom->coordinate(), atomCoordinates(atom->getNeighbors()), 1.0);
                    std::unique_ptr<Atom> newAtom = std::make_unique<Atom>("HHH", newCoord);
                    newAtom->makeInvisible();
                    addBond(atom, newAtom.get());
                    residue->addAtom(std::move(newAtom));
                }
            }
        }
    }
}
