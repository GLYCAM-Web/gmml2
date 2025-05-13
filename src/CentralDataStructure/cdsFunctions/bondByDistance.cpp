#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/logging.hpp"

namespace
{
    struct ResidueAndAtoms
    {
        std::string residueId;
        cds::Residue* residue;
        std::vector<cds::Atom*> atoms;
    };

    std::vector<ResidueAndAtoms> residueAndAtomVector(const pdb::PdbData& pdbData, std::vector<cds::Residue*>& residues)
    {
        std::vector<ResidueAndAtoms> result;
        result.reserve(residues.size());
        for (auto& res : residues)
        {
            size_t index         = codeUtils::indexOf(pdbData.objects.residues, res);
            std::string stringId = pdb::residueStringId(pdbData, index);
            result.push_back({stringId, res, res->getAtoms()});
        }
        return result;
    }

    bool bondAtomsIfClose(cds::Atom* atom1, cds::Atom* atom2)
    {
        double maxLength = MolecularMetadata::maxBondLengthByAtomType(atom1->cachedElement(), atom2->cachedElement());
        if (withinDistance(maxLength, atom1->coordinate(), atom2->coordinate()))
        {
            addBond(atom1, atom2);
            return true;
        }
        return false;
    }

    void bondAtomsAndResiduesByDistance(ResidueAndAtoms& a, ResidueAndAtoms& b)
    {
        cds::Residue* residueA    = a.residue;
        cds::Residue* residueB    = b.residue;
        bool residuesAreConnected = false;
        for (auto& atomA : a.atoms)
        {
            for (auto& atomB : b.atoms)
            {
                if (bondAtomsIfClose(atomA, atomB))
                {
                    residuesAreConnected = true; // only needs to be true once to connect residues.
                }
            }
        }
        if (residuesAreConnected)
        {
            std::string edgeName = a.residueId + "--" + b.residueId;
            residueA->addNeighbor(edgeName, residueB);
        }
    }

    void bondResidueAtoms(std::vector<ResidueAndAtoms>& vec, size_t n)
    {
        std::vector<cds::Atom*> res1Atoms = vec[n].atoms;
        for (size_t k = n + 1; k < vec.size(); k++)
        {
            std::vector<cds::Atom*> res2Atoms = vec[k].atoms;
            if (cds::withinDistance(constants::residueDistanceOverlapCutoff, res1Atoms.at(0)->coordinate(),
                                    res2Atoms.at(0)->coordinate()))
            {
                bondAtomsAndResiduesByDistance(vec[n], vec[k]);
            }
        }
    }
} // namespace

void cds::bondAtomsByDistance(std::vector<cds::Atom*> atoms)
{
    std::vector<MolecularMetadata::Element> elements;
    elements.reserve(atoms.size());
    std::vector<Coordinate> coordinates = atomCoordinates(atoms);
    for (auto& a : atoms)
    {
        elements.push_back(a->cachedElement());
    }
    for (size_t n = 0; n < atoms.size(); n++)
    {
        for (size_t k = n + 1; k < atoms.size(); k++)
        {
            double maxLength = MolecularMetadata::maxBondLengthByAtomType(elements[n], elements[k]);
            if (withinDistance(maxLength, coordinates[n], coordinates[k]))
            {
                addBond(atoms[n], atoms[k]);
            }
        }
    }
}

void cds::bondAtomsAndResiduesByDistance(pdb::PdbData& pdbData, std::vector<cds::Residue*> residues)
{
    std::vector<ResidueAndAtoms> vec = residueAndAtomVector(pdbData, residues);
    for (size_t n = 0; n < residues.size(); n++)
    { // First bond by distance for atoms within each residue
        cds::bondAtomsByDistance(vec[n].atoms);
        // Then for each residue, find other residues within reasonable residue distance.
        bondResidueAtoms(vec, n);
    }
}

void cds::distanceBondIntra(std::vector<cds::Residue*> residues)
{
    for (auto& res : residues)
    {
        cds::bondAtomsByDistance(res->getAtoms());
    }
}

void cds::distanceBondInter(pdb::PdbData& pdbData, std::vector<cds::Residue*> residues)
{
    std::vector<ResidueAndAtoms> vec = residueAndAtomVector(pdbData, residues);
    for (size_t n = 0; n < vec.size(); n++)
    {
        bondResidueAtoms(vec, n);
    }
}
