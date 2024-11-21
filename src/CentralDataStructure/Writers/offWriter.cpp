#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>
#include <string>

namespace
{
    std::vector<size_t> atomIndices(const std::vector<cds::Atom*>& atoms, const std::vector<cds::Atom*>& reindex)
    {
        std::vector<size_t> result;
        result.reserve(reindex.size());
        for (auto& find : reindex)
        {
            result.push_back(cds::atomVectorIndex(atoms, find));
        }
        return result;
    }

    std::vector<std::pair<size_t, size_t>> uniqueAtomBonds(const std::vector<cds::Atom*>& atoms)
    {
        std::vector<std::pair<size_t, size_t>> result;
        result.reserve(2 * atoms.size());
        for (auto& atom : atoms)
        {
            size_t index = cds::atomVectorIndex(atoms, atom);
            for (auto& neighbor : atom->getChildren())
            {
                size_t neighborIndex = cds::atomVectorIndex(atoms, neighbor);
                result.push_back({index, neighborIndex});
            }
        }
        return result;
    }
} // namespace

cds::OffFileData cds::toOffFileData(const std::vector<Residue*>& residues)
{
    std::vector<Atom*> atoms;
    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<size_t>> connections;
    std::vector<size_t> atomResidues;
    size_t residueIndex = 0;
    for (auto& residue : residues)
    {
        std::vector<Atom*> residueAtoms = residue->getAtoms();
        indices.push_back(codeUtils::indexVectorWithOffset(atoms.size(), residueAtoms));
        std::vector<Atom*> atomsConnected = atomsConnectedToOtherResidues(residueAtoms);
        codeUtils::insertInto(atoms, residueAtoms);
        codeUtils::insertInto(atomResidues, std::vector<size_t>(residueAtoms.size(), residueIndex));
        connections.push_back(atomIndices(atoms, atomsConnected));
        residueIndex++;
    }
    OffFileAtomData atomData {atomNumbers(atoms), atomNames(atoms),       atomTypes(atoms), atomAtomicNumbers(atoms),
                              atomCharges(atoms), atomCoordinates(atoms), atomResidues,     uniqueAtomBonds(atoms)};
    OffFileResidueData residueData {residueNumbers(residues), residueNames(residues), residueTypes(residues), indices,
                                    connections};
    return OffFileData {residueData, atomData};
}

void cds::serializeResiduesIndividually(std::vector<cds::Residue*>& residues)
{
    for (auto& residue : residues)
    {
        residue->setNumber(1);
        cds::serializeNumbers(residue->getAtoms());
    }
}
