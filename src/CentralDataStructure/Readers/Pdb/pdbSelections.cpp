#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp"

#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"

using pdb::residueSelector;

std::vector<cds::Atom*> pdb::getAtoms(const std::vector<cds::Assembly*>& assemblies)
{
    std::vector<cds::Atom*> result;
    for (auto& assembly : assemblies)
    {
        codeUtils::insertInto(result, assembly->getAtoms());
    }
    return result;
}

std::vector<cds::Residue*> pdb::getResidues(const std::vector<cds::Assembly*>& assemblies)
{
    std::vector<cds::Residue*> result;
    for (auto& assembly : assemblies)
    {
        codeUtils::insertInto(result, assembly->getResidues());
    }
    return result;
}

size_t pdb::residueSelector(const PdbData& data, const pdb::ResidueId& residueId, const int modelNumber)
{
    for (size_t n = 0; n < data.indices.assemblyCount; n++)
    {
        if (data.objects.assemblies[n]->getNumber() == modelNumber)
        {
            std::vector<size_t> residueAssembly =
                codeUtils::indicesToValues(data.indices.moleculeAssembly, data.indices.residueMolecule);
            return residueSelector(data, codeUtils::indicesOfElement(residueAssembly, n), residueId);
        }
    }
    return data.indices.residueCount;
}

size_t pdb::residueSelector(const PdbData& data, std::vector<size_t> residueIds, const pdb::ResidueId& queryId)
{ // I'm using empty() to mean that it could be anything.
    for (size_t residueId : residueIds)
    {
        ResidueId id = pdbResidueId(data, residueId);
        if (queryId.getName().empty() || queryId.getName() == id.getName())
        {
            if (queryId.getNumber().empty() || queryId.getNumber() == id.getNumber())
            {
                if (queryId.getInsertionCode().empty() || queryId.getInsertionCode() == id.getInsertionCode())
                {
                    if ((queryId.getChainId().empty() || queryId.getChainId() == id.getChainId()))
                    {
                        return residueId;
                    }
                }
            }
        }
    }
    return data.indices.residueCount;
}
