#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"
#include "includes/CodeUtils/logging.hpp"

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
            double maxLength = MolecularMetadata::getMaxBondLengthByAtomType(elements[n], elements[k]);
            if (withinDistance(maxLength, *coordinates[n], *coordinates[k]))
            {
                atoms[n]->addBond(atoms[k]);
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
            if (MolecularMetadata::bondAtomsIfClose(atomA, atomB))
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
        cds::Residue* res1                = *it1;
        std::vector<cds::Atom*> res1Atoms = res1->getAtoms();
        cds::bondAtomsByDistance(res1Atoms);
        // Then for each residue, find other residues within reasonable residue distance.
        for (std::vector<cds::Residue*>::iterator it2 = std::next(it1); it2 != residues.end(); ++it2)
        {
            cds::Residue* res2                = *it2;
            std::vector<cds::Atom*> res2Atoms = res2->getAtoms();
            double residueDistance            = res1Atoms.at(0)->calculateDistance(res2Atoms.at(0));
            if (residueDistance < constants::residueDistanceOverlapCutoff)
            {
                cds::bondAtomsAndResiduesByDistance(res1, res2);
            }
        }
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
        cds::Residue* res1                = *it1;
        std::vector<cds::Atom*> res1Atoms = res1->getAtoms();
        for (std::vector<cds::Residue*>::iterator it2 = std::next(it1); it2 != residues.end(); ++it2)
        {
            cds::Residue* res2                = *it2;
            std::vector<cds::Atom*> res2Atoms = res2->getAtoms();
            double residueDistance            = res1Atoms.at(0)->calculateDistance(res2Atoms.at(0));
            if (residueDistance < constants::residueDistanceOverlapCutoff)
            {
                cds::bondAtomsAndResiduesByDistance(res1, res2);
            }
        }
    }
}
