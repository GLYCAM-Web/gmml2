#include "includes/Graph/cdsToGraph.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace graph
{
    GraphDataLayer toGraphDataLayer(cds::Assembly& assembly)
    {
        GraphDataLayer result;
        std::vector<cds::Residue*> residues;
        std::vector<size_t> residueMolecule;
        std::vector<cds::Atom*> atoms;
        std::vector<size_t> atomResidue;

        for (cds::Molecule* molecule : assembly.getMolecules())
        {
            size_t moleculeId = result.addMolecule(molecule->getName());
            for (cds::Residue* residue : molecule->getResidues())
            {
                size_t residueId = result.addResidue(moleculeId, residue->getName());
                for (cds::Atom* atom : residue->getAtoms())
                {
                    result.addAtom(residueId, atom->coordinate(), atom->getName(), atom->getElement(),
                                   atom->getCharge());
                    atoms.push_back(atom);
                    atomResidue.push_back(residueId);
                }
                residues.push_back(residue);
                residueMolecule.push_back(moleculeId);
            }
        }

        for (size_t n = 0; n < atoms.size(); n++)
        {
            cds::Atom* atom = atoms[n];
            for (cds::Atom* child : atom->getChildren())
            {
                size_t k          = codeUtils::indexOf(atoms, child);
                bool sameResidue  = atomResidue[n] == atomResidue[k];
                bool sameMolecule = residueMolecule[atomResidue[n]] == residueMolecule[atomResidue[k]];
                std::optional<ResidueLinkageStruct> residueLinkage =
                    sameResidue ? std::nullopt
                                : std::optional<ResidueLinkageStruct>(ResidueLinkageStruct {std::vector<size_t>()});
                std::optional<MoleculeLinkageStruct> moleculeLinkage =
                    sameMolecule ? std::nullopt : std::optional<MoleculeLinkageStruct>(MoleculeLinkageStruct {});
                result.addBond(n, k, BondType::covalent, residueLinkage, moleculeLinkage);
            }
        }

        return result;
    }
} // namespace graph
