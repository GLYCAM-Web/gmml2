#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp"

#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"

using pdb::residueSelector;

std::vector<cds::Atom*> pdb::getAtoms(const std::vector<cds::Assembly>& assemblies)
{
    std::vector<cds::Atom*> result;
    for (auto& assembly : assemblies)
    {
        codeUtils::insertInto(result, assembly.getAtoms());
    }
    return result;
}

std::vector<cds::Residue*> pdb::getResidues(const std::vector<cds::Assembly>& assemblies)
{
    std::vector<cds::Residue*> result;
    for (auto& assembly : assemblies)
    {
        codeUtils::insertInto(result, assembly.getResidues());
    }
    return result;
}

pdb::PdbResidue* pdb::residueSelector(const PdbFile& pdbFile, const pdb::ResidueId& residueId, const int modelNumber)
{
    for (auto& model : pdbFile.getAssemblies())
    {
        // std::cout << model->getNumber() << ":" << modelNumber << std::endl;
        if (model.getNumber() == modelNumber)
        {
            return pdb::residueSelector(model.getResidues(), residueId);
        }
    }
    return nullptr;
}

pdb::PdbResidue* pdb::residueSelector(std::vector<cds::Residue*> residues, const pdb::ResidueId& queryId)
{ // I'm using empty() to mean that it could be anything.
    for (auto& residue : residues)
    {
        PdbResidue* pdbResidue = codeUtils::erratic_cast<PdbResidue*>(residue);
        // std::cout << "currentId vs queryId: " << pdbResidue->getId() << " vs " << queryId << std::endl;
        if (queryId.getName().empty() || queryId.getName() == pdbResidue->getId().getName())
        {
            if (queryId.getNumber().empty() || queryId.getNumber() == pdbResidue->getId().getNumber())
            {
                if (queryId.getInsertionCode().empty() ||
                    queryId.getInsertionCode() == pdbResidue->getId().getInsertionCode())
                {
                    if ((queryId.getChainId().empty() || queryId.getChainId() == pdbResidue->getId().getChainId()))
                    {
                        return pdbResidue;
                    }
                }
            }
        }
    }
    return nullptr;
}
