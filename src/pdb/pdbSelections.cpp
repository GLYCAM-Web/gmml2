#include "include/pdb/pdbSelections.hpp"

#include "include/CentralDataStructure/assembly.hpp"
#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbModel.hpp"
#include "include/pdb/pdbResidue.hpp"
#include "include/util/containers.hpp"

namespace gmml
{
    namespace pdb
    {
        std::vector<Atom*> getAtoms(const std::vector<Assembly*>& assemblies)
        {
            std::vector<Atom*> result;
            for (auto& assembly : assemblies)
            {
                util::insertInto(result, assembly->getAtoms());
            }
            return result;
        }

        std::vector<Residue*> getResidues(const std::vector<Assembly*>& assemblies)
        {
            std::vector<Residue*> result;
            for (auto& assembly : assemblies)
            {
                util::insertInto(result, assembly->getResidues());
            }
            return result;
        }

        size_t residueSelector(const PdbData& data, const ResidueId& residueId, const int modelNumber)
        {
            for (size_t n = 0; n < data.indices.assemblyCount; n++)
            {
                if (data.objects.assemblies[n]->getNumber() == modelNumber)
                {
                    return residueSelector(data, assemblyResidues(data.indices, n), residueId);
                }
            }
            return data.indices.residueCount;
        }

        size_t residueSelector(const PdbData& data, std::vector<size_t> residueIds, const ResidueId& queryId)
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
    } // namespace pdb
} // namespace gmml
