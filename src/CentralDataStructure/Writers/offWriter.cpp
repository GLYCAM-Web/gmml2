#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
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

    std::vector<cds::Atom*> atomsConnectedToOtherResidues(std::vector<cds::Atom*> atoms)
    {
        std::vector<cds::Atom*> foundAtoms;
        for (auto& atom : atoms)
        {
            bool connected = false;
            for (auto& neighbor : atom->getNeighbors())
            { // check if neighbor is not one of the atoms in this residue.
                connected = connected || !codeUtils::contains(atoms, neighbor);
            }
            if (connected)
            {
                foundAtoms.push_back(atom);
            }
        }
        return foundAtoms;
    }
} // namespace

cds::OffFileData cds::toOffFileData(const std::vector<Residue*>& residues)
{
    std::vector<Atom*> atoms;
    std::vector<std::vector<size_t>> connections;
    for (auto& residue : residues)
    {
        std::vector<Atom*> residueAtoms   = residue->getAtoms();
        std::vector<Atom*> atomsConnected = atomsConnectedToOtherResidues(residueAtoms);
        codeUtils::insertInto(atoms, residueAtoms);
        connections.push_back(atomIndices(atoms, atomsConnected));
    }
    OffFileAtomData atomData {atomNumbers(atoms),       atomNames(atoms),   atomTypes(atoms),
                              atomAtomicNumbers(atoms), atomCharges(atoms), atomCoordinates(atoms)};
    OffFileResidueData residueData {residueNumbers(residues), residueNames(residues), residueTypes(residues),
                                    connections};
    OffFileFormat format;
    return OffFileData {format, residueData, atomData};
}

void cds::serializeResiduesIndividually(std::vector<cds::Residue*>& residues)
{
    for (auto& residue : residues)
    {
        residue->setNumber(1);
        cds::serializeNumbers(residue->getAtoms());
    }
}

void cds::WriteOff(std::ostream& stream, const std::string& name, const GraphIndexData& data)
{
    assembly::Graph graph    = createVisibleAssemblyGraph(data);
    cds::OffFileData offData = cds::toOffFileData(data.objects.residues);
    offData.atoms.numbers    = serializedNumberVector(graph.atoms.source.nodeAlive);
    offData.residues.numbers = serializedNumberVector(data.objects.residues.size());
    cds::WriteResiduesTogetherToOffFile(stream, graph, offData, codeUtils::indexVector(data.objects.residues), name);
}
