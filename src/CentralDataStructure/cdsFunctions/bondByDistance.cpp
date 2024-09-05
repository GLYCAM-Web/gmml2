#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CodeUtils/logging.hpp"

namespace
{
    bool bondAtomsIfClose(cds::Atom* atom1, cds::Atom* atom2)
    {
        double maxLength = MolecularMetadata::maxBondLengthByAtomType(
            MolecularMetadata::toElement(atom1->getElement()), MolecularMetadata::toElement(atom2->getElement()));
        if (withinDistance(maxLength, *atom1->getCoordinate(), *atom2->getCoordinate()))
        {
            addBond(atom1, atom2);
            return true;
        }
        return false;
    }

    void bondResidueAtoms(std::vector<cds::Residue*>::iterator it1, std::vector<cds::Residue*>::iterator itEnd)
    {
        cds::Residue* res1                = *it1;
        std::vector<cds::Atom*> res1Atoms = res1->getAtoms();
        for (std::vector<cds::Residue*>::iterator it2 = std::next(it1); it2 != itEnd; ++it2)
        {
            cds::Residue* res2                = *it2;
            std::vector<cds::Atom*> res2Atoms = res2->getAtoms();
            double residueDistance = distance(*res1Atoms.at(0)->getCoordinate(), *res2Atoms.at(0)->getCoordinate());
            if (residueDistance < constants::residueDistanceOverlapCutoff)
            {
                cds::bondAtomsAndResiduesByDistance(res1, res2);
            }
        }
    }
} // namespace

void cds::bondAtomsByDistance(std::vector<cds::Atom*> atoms)
{
    std::vector<MolecularMetadata::Element> elements;
    elements.reserve(atoms.size());
    std::vector<Coordinate*> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& a : atoms)
    {
        elements.push_back(MolecularMetadata::toElement(a->getElement()));
        coordinates.push_back(a->getCoordinate());
    }
    for (size_t n = 0; n < atoms.size(); n++)
    {
        for (size_t k = n + 1; k < atoms.size(); k++)
        {
            double maxLength = MolecularMetadata::maxBondLengthByAtomType(elements[n], elements[k]);
            if (withinDistance(maxLength, *coordinates[n], *coordinates[k]))
            {
                addBond(atoms[n], atoms[k]);
            }
        }
    }
}

void cds::bondAtomsAndResiduesByDistance(cds::Residue* residueA, cds::Residue* residueB)
{
    bool residuesAreConnected = false;
    for (auto& atomA : residueA->getAtoms())
    {
        for (auto& atomB : residueB->getAtoms())
        {
            if (bondAtomsIfClose(atomA, atomB))
            {
                residuesAreConnected = true; // only needs to be true once to connect residues.
            }
        }
    }
    if (residuesAreConnected)
    {
        std::string edgeName = residueA->getStringId() + "--" + residueB->getStringId();
        residueA->addNeighbor(edgeName, residueB);
    }
}

void cds::bondAtomsAndResiduesByDistance(std::vector<cds::Residue*> residues)
{
    for (std::vector<cds::Residue*>::iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
    { // First bond by distance for atoms within each residue
        cds::bondAtomsByDistance((*it1)->getAtoms());
        // Then for each residue, find other residues within reasonable residue distance.
        bondResidueAtoms(it1, residues.end());
    }
}

void cds::distanceBondIntra(std::vector<cds::Residue*> residues)
{
    for (auto& res : residues)
    {
        cds::bondAtomsByDistance(res->getAtoms());
    }
}

void cds::distanceBondInter(std::vector<cds::Residue*> residues)
{
    for (std::vector<cds::Residue*>::iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
    {
        bondResidueAtoms(it1, residues.end());
    }
}
