#include "include/CentralDataStructure/residueLinkage/psiAngleHydrogen.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/measurements.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <memory>
#include <vector>

namespace gmml
{
    Atom* findHydrogenForPsiAngle(const Atom* atom)
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

    void createHydrogenForPsiAngles(
        const DihedralAngleDataTable& metadataTable,
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
} // namespace gmml
